
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

from lsst.daf.persistence import Butler

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from scipy.spatial import cKDTree


class MakeDecapsMatchedCatalog:


    def run(self):
        parser = argparse.ArgumentParser()
        parser.add_argument("visit", type=int)
        parser.add_argument("repo", type=str)
        parser.add_argument("crowdsource_catalog", type=str)

        args = parser.parse_args()

        butler = Butler(args.repo)
        self.matchVisit(butler, args.visit,
                        crowdsource_filename=args.crowdsource_catalog)


    def matchVisit(self, butler, visitId, crowdsource_filename):

        # crowdsource_filename = "catalogs_decaps/c4d_160317_001229_ooi_g_v1.cat.fits"
        camera = butler.get("camera")

        match_tables = []
        for n in range(1, 63):
            if(not butler.datasetExists("crowdedsrc", visit=visitId, ccd=n)):
                continue
            dataRef = butler.dataRef("crowdedsrc", visit=visitId, ccd=n)
            sensor_match = self.matchDetector(dataRef, crowdsource_filename)
            match_tables.append(sensor_match)

        combined_table = pd.concat(match_tables)
        combined_table.to_parquet("combined_visit_{:d}.parquet".format(visitId))


    def matchDetector(self, sensorRef, crowdsource_filename, max_dist_pixels=1.5):
        detector = sensorRef.get("calexp_detector")

        detector_name = detector.getName()

        crowdsource_file = fits.open(crowdsource_filename)

        crowdsource_table = Table.read(crowdsource_filename,
                                       hdu=f"{detector_name}_CAT")

        crowd_tree = cKDTree(np.stack((crowdsource_table['x'],
                                       crowdsource_table['y']), axis=1))

        stack_table = sensorRef.get("crowdedsrc")
        stack_tree = cKDTree(np.stack((stack_table['centroid_y'],
                                       stack_table['centroid_x']), axis=1))

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


