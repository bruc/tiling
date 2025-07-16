README  (Last Updated 09/02/2024)

********************
INTRODUCTION
********************

This README file describes the contents in this directory.

This dataset contains raw, reprocessed, and aligned files for multiple AAV samples.

The 2021 and 2022 libraries were sequenced on the Sequel® II System and processed using SMRT® Link v10.1 or SMRT Link v11.0 [1], followed by community tool analysis [2].

The 2023-Revio-ssAAV-pAV-CMF-GFP library was sequenced on the Revio system and processed using a preliminary version of SMRT Link v13.0, followed by FormBio’s LAAVA pipeline (version 2.5.0) .

To learn more about PacBio's gene therapy research tools offering, see: https://www.pacb.com/gene-therapy/


********************
SAMPLE
********************

2021-scAAV/rh10-CBA-eGFP 
Vendor – Vector Biolabs
Lot No – 210621-210723

2022-pAV-CMV-GFP
Serotype - AAV9
Vendor – Vigene Biosciences
Catalog No – CV10009 
Lot No (Date) - 2021.02.25

2023-Revio-ssAAV-pAV-CMF-GFP
Serotype - AAV9
Vendor - Charles River
Catalog No - CV10009 Control GFP AAV9 Virus
Lot No - 3A218-C-1

********************
METHODS
********************
DNA extraction – PureLink Viral DNA/RNA Mini Kit

Library Preparation: 

Procedure & Checklist - Preparing multiplexed AAV SMRTbell® libraries using SMRTbell prep kit 3.0

Sequencing: 

2021-scAAV-CBA-eGFP
Sequel II system with Sequel II binding kit 2.1 and Sequel II sequencing kit 2.0 (4 rxn), Sequencing primer v4

2022-ssAAV-pAV-CMF-GFP
Sequel II system with Sequel II binding kit 2.0 and Sequel II sequencing kit 2.0 (4 rxn), Sequencing primer v4

2022-ssAAV-scAAV-mix
1 part of ssAAV library mixed with 10 parts of scAAV library
Sequel II system with Sequel II binding kit 3.1 and Sequel II sequencing kit 2.0 (4 rxn), Sequencing primer v4

2023-Revio-ssAAV-pAV-CMF-GFP
Revio system with Revio polymerase kit and Revio sequencing plate. 


Analysis: 

Refer to the LAAVA GitHub (https://github.com/formbio/LAAVA)
and Talevich et al. 2024 (https://www.biorxiv.org/content/10.1101/2024.05.07.592296v1)

   
********************
FILE DESCRIPTION
********************

Each sample will contain the following folders:

========================
0-reads
========================

CCS reads after running "AAV" mode in Run Design via SMRT Link

0-reads/
|---- <movie>.hifi_reads.bam
|---- <movie>.hifi_reads.bam.pbi



========================
1-(demuxed)-and-mapped
========================

If samples are multiplexed, they are shown as demux individual samples here.
For each sample, the HiFi reads are mapped to the reference plasmid.
The key files are shown below.

1-mapped/
|---- <sample>_AAV_report.pdf
|---- <sample>.nonmatch_stat.csv.gz
|---- <sample>.per_read.csv
|---- <sample>.summary.csv
|---- <sample>.tagged.sorted.bam
|---- <sample>.tagged.sorted.bam.bai


This aligned sorted BAM file is what you can visualize using IGV.




4. REFERENCES

[1] Procedure & checklist: Preparing multiplexed AAV SMRTbell® libraries using SMRTbell prep kit 3.0
https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-multiplexed-AAV-SMRTbell-libraries-using-SMRTbell-prep-kit-3.0.pdf

[2] AAV GitHub Wiki: https://github.com/formbio/laava


For Research Use Only. Not for use in diagnostic procedures.  Copyright 2024, 
Pacific Biosciences of California, Inc. All rights reserved. The data provided in 
these files is subject to change without notice and Pacific Biosciences assumes no 
responsibility for any errors or omissions. Certain notices, terms, conditions and/or 
use restrictions may pertain to your use of Pacific Biosciences data, products and/or 
third party products. Please refer to the applicable Pacific Biosciences Terms and 
Conditions of Sale and to the applicable license terms at 
http://www.pacificbiosciences.com/licenses.html.

