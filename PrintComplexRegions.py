#!/usr/bin/env python

import argparse
import sys
ap = argparse.ArgumentParser(description="Print regions from gap bed files that are complex.")
ap.add_argument("bed", help="Input bed file.  The annotation is the 7th field.")
ap.add_argument("out", help="Output file.")
ap.add_argument("-s", "--simple", help="Print simple events instead of complex.", action='store_true', default=False)
ap.add_argument("-1", "--single", help="Print single .", action='store_true', default=False)
ap.add_argument("-i", "--intermediate", help="Print intermediate help files instead of simple or complex, such as double-Alu insertions.", action='store_true', default=False)
ap.add_argument("-a", "--annotation", help="Use this column for repeat annotation", type=int, default=6)

args = ap.parse_args()

action = 'all'

if (args.simple):
    action = 'simple'
if (args.intermediate):
    action = 'intermediate'
if (args.single):
    action = 'single'

gapBed = open(args.bed)
outGapBed = open(args.out, 'w')
reps = ["L2", "5S","7SK","7SLRNA","7SLRNA","ACRO1","ACRO1","ALR","ALRa","ALRa","ALRb","ALRb","ALR","AluJb","AluJo","AluS", "AluSc","AluSg1","AluSg","AluSp","AluSq","AluSx","AluSx","AluYa5","AluYa8","AluYa8","AluYb8","AluYb9","AluYc3","AluYc5","AluYc","AluYd8","AluYg6","AluYh9","AluY","AluY","Arthur1","BC200","BLACKJACK","BSRa","BSRb","BSRd","BSRf","BSR","CER","Charlie10","Charlie11","Charlie12","Charlie1a","Charlie1b","Charlie1","Charlie2a","Charlie2b","Charlie3","Charlie4a","Charlie4","Charlie5","Charlie6","Charlie7","Charlie8","Charlie9","Cheshire","D20S16","D20S16","FAM","FAM","FLAM", "FLAM_C", "FordPrefect","FordPrefect","FRAM","GSAT","GSATII","GSATII","GSAT","GSATX","GSATX","HAL1b","HAL1b","HAL1","HAL1","HSAT4","HSAT5","HSAT5","HSAT6","HSATI","HSMAR1","HSMAR2","HY1","HY3","HY4","HY5","Kanga1","Kanga2","L1","L1HS","L1","L1M1","L1M2","L1M2a","L1M2a1","L1M2b","L1M2c", "L1M3", "L1M3a","L1M3b","L1M3c","L1M3d","L1M3de","L1M3e","L1M3f","L1M4","L1M4b","L1M4c","L1MA1","L1MA10","L1MA2","L1MA3","L1MA4","L1MA4A","L1MA5","L1MA5A","L1MA6","L1MA7","L1MA8","L1MA9","L1MB1","L1MB2","L1MB3","L1MB4","L1MB5","L1MB7","L1MB8","L1MC1","L1MC2","L1MC3","L1MC4","L1MC4a","L1MC5","L1MCa","L1MCb","L1MCc","L1MD1","L1MD2","L1MD3","L1MDa","L1MDb","L1ME1","L1ME2","L1ME3","L1ME3A","L1ME3B","L1ME4a","L1MEa","L1MEb","L1MEc","L1MEd","L1MEe","L1P3","L1P3b","L1PA10","L1PA11","L1PA12","L1PA13","L1PA13","L1PA14","L1PA15","L1PA15-16","L1PA16","L1PA16","L1PA17","L1PA2","L1PA3","L1PA4","L1PA4","L1PA5","L1PA6","L1PA7","L1PA7","L1PA8","L1PA8A","L1PB1","L1PB1","L1PB2","L1PB3","L1PB4","L1PB4","L1PBa","L1PBa1","L1PBb", "L1P", "Looper","LOR1a","LOR1b","LSAU","LTR10A","LTR10A","LTR10B1","LTR10B","LTR10C","LTR10D","LTR10E","LTR10F","LTR10G","LTR11","LTR13A","LTR13","LTR14A","LTR14B","LTR14C","LTR14","LTR15","LTR16A1","LTR16A","LTR16B","LTR16C","LTR16D","LTR18A","LTR18B","LTR19A","LTR19B","LTR19C","LTR1B","LTR1C","LTR1D","LTR1","LTR21A","LTR21B","LTR22A","LTR22B","LTR22C","LTR22","LTR23","LTR24B","LTR24C","LTR24","LTR25","LTR26B","LTR26E","LTR26","LTR27B","LTR27","LTR28","LTR29","LTR2B","LTR2C","LTR2","LTR2","LTR30","LTR31","LTR32","LTR33A","LTR33","LTR34","LTR35","LTR36","LTR37A","LTR37B","LTR38B","LTR38C","LTR38","LTR39","LTR3A","LTR3B","LTR3B","LTR3","LTR40a","LTR40b","LTR40c","LTR41","LTR42","LTR43B","LTR43","LTR44","LTR45B","LTR45C","LTR45","LTR46","LTR47A","LTR47B","LTR48B","LTR48","LTR49","LTR4","LTR5","LTR50","LTR51","LTR52","LTR53","LTR54B","LTR54","LTR55","LTR56","LTR57","LTR58","LTR59","LTR5A","LTR5B","LTR5B","LTR5","LTR5","LTR60","LTR61","LTR62","LTR64","LTR65","LTR66","LTR67","LTR68","LTR69","LTR6A","LTR6B","LTR70","LTR71A","LTR71B","LTR72B","LTR72","LTR73","LTR75","LTR75","LTR76","LTR77","LTR7A","LTR7B","LTR7","LTR7","LTR8A","LTR8","LTR9B","LTR9","MADE1","MADE2","MARNA","MER101B","MER101","MER102a","MER102b","MER103","MER104A","MER104B","MER104C","MER104","MER105","MER106A","MER106B","MER107","MER109","MER110A","MER110","MER112","MER113","MER115","MER117","MER119","MER11A","MER11B","MER11C","MER11C","MER11D","MER121","MER1A","MER1B","MER20B","MER20","MER21A","MER21B","MER21C","MER2B","MER2","MER30B","MER30","MER31A","MER31B","MER33","MER34B","MER34C","MER34C","MER34D","MER34","MER39B","MER39","MER3","MER41A","MER41A","MER41B","MER41C","MER41D","MER41E","MER41G","MER44A","MER44B","MER44C","MER44D","MER45A","MER45B","MER45C","MER45R","MER46C","MER47A","MER47B","MER47C","MER48","MER49","MER4A1","MER4A","MER4A","MER4B","MER4C","MER4D0","MER4D1","MER4D","MER4E1","MER4E","MER50B","MER50","MER51A","MER51A","MER51B","MER51C","MER51D","MER51E","MER52A","MER52C","MER52D","MER53","MER54A","MER54B","MER57A","MER57B","MER58A","MER58B","MER58C","MER58D","MER5A1","MER5A","MER5B","MER5C","MER61A","MER61B","MER61C","MER61D","MER61E","MER61F","MER63A","MER63B","MER63C","MER63D","MER65A","MER65B","MER65C","MER65D","MER66A","MER66B","MER66C","MER67A","MER67B","MER67C","MER67D","MER68","MER69A","MER69B","MER6A","MER6B","MER6C","MER6","MER70A","MER70B","MER72B","MER72","MER73","MER74A","MER74B","MER74C","MER75B","MER75","MER76","MER77","MER81","MER82","MER83B","MER83C","MER83","MER84","MER85","MER85","MER87","MER88","MER89","MER8","MER90a","MER90","MER91A","MER91B","MER91C","MER92A","MER92B","MER92C","MER93a","MER93b","MER93B","MER94","MER95","MER96B","MER96","MER97a","MER97b","MER97c","MER99","MER9B","MER9","Merlin1","MIR3","MIR","MLT1A0","MLT1A1","MLT1A","MLT1B","MLT1C","MLT1D","MLT1E1A","MLT1E1","MLT1E2","MLT1E3","MLT1E","MLT1F1","MLT1F2","MLT1F-int","MLT1F","MLT1G1","MLT1G3","MLT1G","MLT1H1","MLT1H2","MLT1H-int","MLT1H","MLT1I","MLT1-int","MLT1J1","MLT1J2","MLT1J","MLT1K","MLT1L","MLT2A1","MLT2A2","MLT2B1","MLT2B2","MLT2B3","MLT2B4","MLT2B5","MLT2C1","MLT2C2","MLT2D","MLT2E","MLT2F","MLT-int","MSR1","MSTA","MSTB1","MSTB2","MSTB","MSTC","MSTD","MST-int","ORSL","PABL","PMER1","PRIMA4","REP522","Ricksha","SAR","SATR1","SATR2","SST1","SUBTEL","SVA","TAR1","THE1A","THE1B","THE1B","THE1C","THE1D","THE1-int","Tigger1","Tigger2a","Tigger2","Tigger3a","Tigger3b","Tigger3c","Tigger3d","Tigger3","Tigger4a","Tigger4b","Tigger4","Tigger5","Tigger6a","Tigger6b","Tigger7","Tigger8","tRNA-Ala-GCA","tRNA-Ala-GCG","tRNA-Ala-GCY","tRNA-Ala-GCY","tRNA-Arg-AGA","tRNA-Arg-AGG","tRNA-Arg-CGA","tRNA-Arg-CGA","tRNA-Arg-CGG","tRNA-Arg-CGY","tRNA-Arg-CGY","tRNA-Asn-AAC","tRNA-Asn-AAT","tRNA-Asp-GAY","tRNA-Cys-TGY","tRNA-Gln-CAA","tRNA-Gln-CAA","tRNA-Gln-CAG","tRNA-Glu-GAA","tRNA-Glu-GAG","tRNA-Glu-GAG","tRNA-Gly-GGA","tRNA-Gly-GGG","tRNA-Gly-GGY","tRNA-His-CAY","tRNA-His-CAY","tRNA-Ile-ATA","tRNA-Ile-ATC","tRNA-Ile-ATT","tRNA-Leu-CTA","tRNA-Leu-CTA","tRNA-Leu-CTG","tRNA-Leu-CTY","tRNA-Leu-TTA(m)","tRNA-Leu-TTA","tRNA-Leu-TTG","tRNA-Lys-AAA","tRNA-Lys-AAG","tRNA-Met","tRNA-Met-i","tRNA-Met","tRNA-Phe-TTY","tRNA-Pro-CCA","tRNA-Pro-CCG","tRNA-Pro-CCY","tRNA-SeC(e)-TGA","tRNA-Ser-AGY","tRNA-Ser-TCA","tRNA-Ser-TCA(m)","tRNA-Ser-TCA","tRNA-Ser-TCG","tRNA-Ser-TCY","tRNA-Thr-ACA","tRNA-Thr-ACG","tRNA-Thr-ACG","tRNA-Thr-ACY","tRNA-Thr-ACY","tRNA-Trp-TGG","tRNA-Tyr-TAC","tRNA-Tyr-TAT","tRNA-Val-GTA","tRNA-Val-GTG","tRNA-Val-GTY","U13","U13","U14","U17","U1","U2","U3","U4","U5","U6","U7","U8","Zaphod2","Zaphod"]


def CountAll(annot):

    base = [a.split(":")[0] for a in annot]
    vals = [Count(annot, rep) for  rep in reps]
    return sum(vals)

def Match(i,j):
    if (i == j):
        return 1
    else:
        return 0
    
def Count(annot, rep):
    base = [a.split(":")[0] for a in annot]
    return sum([Match(a, rep) for a in base])

for line in gapBed:
    vals = line.split()
    annotations = vals[args.annotation].split(',')
    nRep = CountAll(annotations)
    
    if (action == 'single'):
        if (nRep <= 1):
            outGapBed.write(line)
    elif (action == 'intermediate'):
        if (len(annotations) < 3):
            outGapBed.write(line)
    else:
        if (args.simple):
            if (nRep <= 2):
                outGapBed.write(line)
        else:
            if (nRep >= 3):
                outGapBed.write(line)

outGapBed.close()        
    

