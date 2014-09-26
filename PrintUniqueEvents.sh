#!/usr/bin/env bash
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $1  --prefix AluY --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder $2/1.bed > $2/AluY.simple.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/1.bed --prefix AluS --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder $2/2.bed > $2/AluS.simple.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/2.bed --prefix NONE --minPrefix 1 --maxPrefix 1 --maxNotPrefix 0 --remainder $2/3.bed > $2/NONE.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/3.bed --minSTR 1  --maxNotPrefix 0 --remainder $2/4.bed > $2/STR.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/4.bed --prefix L1HS  --maxNotPrefix 0 --remainder $2/5.bed > $2/L1HS.simple.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/5.bed --prefix Alu  --minPrefix 1 --maxNotPrefix 0 --maxSTR 0 --remainder $2/6.bed > $2/Alu.Mosaic.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/6.bed --prefix Alu  --minSTR 1 --minPrefix 1 --maxNotPrefix 0 --remainder $2/7.bed > $2/Alu.STR.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/7.bed --prefix ALR   --minPrefix 1 --maxNotPrefix 0 --remainder $2/8.bed > $2/ALR.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/8.bed --prefix SVA   --minPrefix 1 --maxNotPrefix 0 --remainder $2/9.bed > $2/SVA.simple.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/9.bed --prefix HERV   --minPrefix 1 --maxNotPrefix 0 --remainder $2/10.bed > $2/HERV.simple.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/10.bed --prefix L1P   --minPrefix 1 --maxNotPrefix 0 --remainder $2/11.bed > $2/L1P.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/11.bed --prefix BSR/Beta   --minPrefix 1 --maxNotPrefix 0 --remainder $2/12.bed > $2/Beta.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/12.bed --prefix HSAT   --minPrefix 1 --maxNotPrefix 0 --remainder $2/13.bed > $2/HSAT.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/13.bed --prefix MER   --minPrefix 1 --maxNotPrefix 0 --remainder $2/14.bed > $2/MER.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/14.bed --prefix L1   --minPrefix 1 --maxNotPrefix 0 --remainder $2/15.bed > $2/L1.bed
/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/14.bed --prefix LTR  --minPrefix 1 --maxNotPrefix 0 --remainder $2/16.bed > $2/LTR.bed$

/net/eichler/vol5/home/mchaisso/projects/PacBioSequencing/scripts/PrintUniqueEvents.py $2/16.bed --max 1 --remainder $2/17.bed > $2/Singletons.bed
mv $2/17.bed $2/Complex.bed


