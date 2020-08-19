
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import pickle
import glob
from astropy.io import ascii

from lsst.daf.persistence import Butler

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree


class MakeDecapsMatchedCatalog:


    def run(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("--use-src-catalog", action="store_true")
        parser.add_argument("visit", type=int)
        parser.add_argument("repo", type=str)
        parser.add_argument("crowdsource_catalog", type=str)

        args = parser.parse_args()

        butler = Butler(args.repo)
        if(args.use_src_catalog):
            self.matchVisit(butler, args.visit,
                            crowdsource_filename=args.crowdsource_catalog,
                            catalog_datatype="src",
                            base_output_name="combined_visit_src")
        else:
            self.matchVisit(butler, args.visit,
                            crowdsource_filename=args.crowdsource_catalog,
                            catalog_datatype="crowdedsrc",
                            base_output_name="combined_visit")


    def matchVisit(self, butler, visitId, crowdsource_filename,
                   base_output_name="combined_visit",
                   catalog_datatype=None):

        # crowdsource_filename = "catalogs_decaps/c4d_160317_001229_ooi_g_v1.cat.fits"
        camera = butler.get("camera")

        match_tables = []
        for n in range(1, 63):
            if(not butler.datasetExists(catalog_datatype, visit=visitId, ccd=n)):
                continue
            dataRef = butler.dataRef(catalog_datatype, visit=visitId, ccd=n)
            sensor_match = self.matchDetector(dataRef, crowdsource_filename,
                                              catalog_datatype=catalog_datatype)
            match_tables.append(sensor_match)

        combined_table = pd.concat(match_tables)
        combined_table.to_parquet("{:s}_{:d}.parquet".format(base_output_name, visitId))


    def matchDetector(self, sensorRef, crowdsource_filename,
                      max_dist_pixels=1.5, catalog_datatype="crowdedsrc"):
        detector = sensorRef.get("calexp_detector")

        detector_name = detector.getName()

        crowdsource_file = fits.open(crowdsource_filename)

        crowdsource_table = Table.read(crowdsource_filename,
                                       hdu=f"{detector_name}_CAT")

        crowd_tree = cKDTree(np.stack((crowdsource_table['x'],
                                       crowdsource_table['y']), axis=1))

        stack_table = sensorRef.get(catalog_datatype)
        stack_tree = cKDTree(np.stack((stack_table['slot_Centroid_y'],
                                       stack_table['slot_Centroid_x']), axis=1))

        matches = stack_tree.query_ball_tree(crowd_tree, max_dist_pixels)
        idx_pairs = list((x,y[0]) for x,y in zip(range(len(stack_table)),
                                                 matches) if len(y) > 0)
        stack_idx = np.array(list(x[0] for x in idx_pairs))
        crowd_idx = np.array(list(x[1] for x in idx_pairs))
        print("Matches: {:d}".format(len(idx_pairs)))

        crowd_unmatched_idx = (set(range(len(crowdsource_table)))
                                - set(x[1] for x in idx_pairs))
        stack_unmatched_idx = (set(range(len(stack_table)))
                                - set(x[0] for x in idx_pairs))

        stack_pandas = stack_table.asAstropy().to_pandas()
        crowd_pandas = crowdsource_table.to_pandas()

        stack_pandas['visit'] = sensorRef.dataId['visit']
        stack_pandas['ccd'] = sensorRef.dataId['ccd']
        stack_pandas['key'] = np.nan
        stack_pandas.loc[stack_idx, 'key'] = crowd_idx
        combined_total = pd.merge(stack_pandas, crowd_pandas, left_on="key",
                                  right_index=True, how="outer")

        return combined_total


class MakeDecapsPlots:

    def run(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("visit", type=int)
        parser.add_argument("catalog", type=str)

        args = parser.parse_args()

        catalog = pd.read_parquet(args.catalog)
        plot_funcs = [self.plot_completeness]
        for f in plot_funcs:
            f(catalog, args.visit)

    def plot_completeness(self, catalog, visit):
        sel_crowdsource = catalog['decapsid'] > 0
        sel_stack = catalog['id'] > 0
        sel_matches =  (catalog['id'] > 0) & (catalog['decapsid'] > 0)
        sel_missing_from_stack = (~(catalog['id'] > 0)) & (catalog['decapsid'] > 0)
        number_ccds = len(np.unique(catalog['ccd'][catalog['ccd'] >= 0].astype(int)))


        H_total, bins = np.histogram(-2.5*np.log10(catalog['flux'][sel_crowdsource]), range=(-18, -6), bins=20)
        H_matches, _ = np.histogram(-2.5*np.log10(catalog['flux'][sel_matches]), bins=bins)
        H_nonmatches, _ = np.histogram(-2.5*np.log10(catalog['flux'][sel_missing_from_stack]), bins=bins)

        ax = plt.gca()

        ax2 = plt.twinx()
        ax2.plot(bins[1:], H_total, 'k--', drawstyle="steps-pre", label="Counts")
        ax2.set_ylabel("Number of stars per bin")

        ax.plot(bins[1:], H_nonmatches/H_total, '-', lw=4,
                drawstyle="steps-pre", label="Fraction missing")
        ax.set_ylabel("Fraction of stars missing from pipe_crowd")
        ax.set_xlabel("Instrumental Mag")

        plt.legend(loc=0, frameon=False)
        plt.savefig(f"matches_visit_{visit}.pdf")

        metrics = {"visit": visit,
                   "processed_ccds": number_ccds,
                   "n_matches": np.sum(sel_matches),
                   "n_stack": np.sum(sel_stack),
                   "n_decaps": np.sum(sel_crowdsource),
                   "n_missing_from_stack": np.sum(sel_missing_from_stack)}

        with open(f"metrics_visit_{visit}.pkl", "wb") as f:
            pickle.dump(metrics, f)


class MakeSummaryPlots:

    def run(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("metricglob", type=str)
        parser.add_argument("--alt-metric-glob", type=str)

        args = parser.parse_args()

        metrics = []
        for filename in glob.glob(args.metricglob):
            with open(filename, "rb") as f:
                metric = pickle.load(f)
                metrics.append(metric)

        alt_metrics = []
        if(len(args.alt_metric_glob) > 0):
            for filename in glob.glob(args.alt_metric_glob):
                with open(filename, "rb") as f:
                    metric = pickle.load(f)
                    alt_metrics.append(metric)

        plot_funcs = [self.plot_all_visit_completeness,
                      self.plot_all_visit_comparitive_counts]
        for f in plot_funcs:
            f(metrics, alt_metrics=alt_metrics)

    def plot_all_visit_comparitive_counts(self, metrics, alt_metrics=None):


        decaps_table = ascii.read("decaps.csv")
        decaps_table_visit = decaps_table['col1']
        decaps_table_density = decaps_table['col6']
        decaps_density_dict = {visit: density for (visit, density) in
                               zip(decaps_table_visit, decaps_table_density)}

        plt.clf()
        ax = plt.gca()

        stack_counts = np.array([x['n_stack']
                              for x in alt_metrics])
        n_ccds = np.array([x['processed_ccds']
                              for x in alt_metrics])
        area = n_ccds * 0.045 # sq degrees per sensor
        decaps_density = np.array([decaps_density_dict[x['visit']] for x in
                                   alt_metrics])
        # ax.plot(stack_counts/area, decaps_counts/area, 'rx', label="stock stack")
        ax.plot(stack_counts/area, decaps_density, 'rx', label="stock stack")

        pipe_crowd_counts = np.array([x['n_stack']
                              for x in metrics])
        n_ccds = np.array([x['processed_ccds']
                              for x in metrics])
        decaps_density = np.array([decaps_density_dict[x['visit']] for x in
                                   metrics])
        area = n_ccds * 0.045 # sq degrees per sensor
        ax.plot(pipe_crowd_counts/area, decaps_density, 'bo', label="pipe_crowd")

        ax.plot([0, 550e3], [0, 550e3], 'k-')

        metrics_by_visit = {x['visit']: x for x in metrics}
        alt_metrics_by_visit = {x['visit']: x for x in alt_metrics}
        for visit in metrics_by_visit.keys():
            original_x = (alt_metrics_by_visit[visit]['n_stack'] /
                          (alt_metrics_by_visit[visit]['processed_ccds'] * 0.045))

            original_y = (alt_metrics_by_visit[visit]['n_decaps'] /
                          (alt_metrics_by_visit[visit]['processed_ccds'] * 0.045))

            new_x = (metrics_by_visit[visit]['n_stack'] /
                     (metrics_by_visit[visit]['processed_ccds'] * 0.045))

            new_y = (metrics_by_visit[visit]['n_decaps'] /
                     (metrics_by_visit[visit]['processed_ccds'] * 0.045))

            # plt.plot([original_x, new_x], [original_y, new_y], 'k-')
            # plt.text(new_x, new_y, "{:d}".format(visit))

        ax.set_xlabel("LSST Density (per sq deg)")
        ax.set_ylabel("DECAPS Density (per sq deg)")
        ax.set_xlim(0, 600e3)
        ax.set_ylim(0, 600e3)
        plt.legend(loc=0, frameon=False)
        plt.savefig("all_visits_comparitive_counts.pdf")


    def plot_all_visit_completeness(self, metrics, alt_metrics=None):
        plt.clf()
        ax = plt.gca()

        fractions = np.array([x['n_matches']/x['n_decaps']
                              for x in metrics])
        decaps_counts = np.array([x['n_decaps']
                              for x in metrics])

        ax.plot(decaps_counts, fractions, 'o')

        if(alt_metrics):
            alt_fractions = np.array([x['n_matches']/x['n_decaps']
                                  for x in alt_metrics])
            alt_decaps_counts = np.array([x['n_decaps']
                                  for x in alt_metrics])
            ax.plot(alt_decaps_counts, alt_fractions, 'ro')

        ax.set_ylim(0, 1.1)
        plt.savefig("all_visit_completness.pdf")




