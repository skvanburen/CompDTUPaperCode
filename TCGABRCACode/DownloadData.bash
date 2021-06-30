#!/bin/bash

#  DownloadData.bash
#  
#
#  Created by Scott Van Buren on 9/20/20.
#  
#SBATCH --time=95:59:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=29g
#SBATCH --output=/nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/DownloadData.out

## Modify the manifest used as needed
## final_manifest_high_purity.txt
## final_manifest_non_high_purity.txt
module load python/3.6.6
/nas/longleaf/home/skvanbur/bin/gdc-client download -m /nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/FilesForGDCDownloader/final_manifest_missing.txt -t /nas/longleaf/home/skvanbur/res/TCGABRCAAnalysis/FilesForGDCDownloader/gdc-user-token.2020-09-20T19_44_18.259Z.txt -d /pine/scr/s/k/skvanbur/TCGABRCAAnalysis/Data/
