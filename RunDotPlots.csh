while read line ; do 
		read=`echo $line | cut -f 1 -d' '`;
		r=`echo $line | cut -f 2 -d' '`;
    a=`echo $line | cut -f 1 -d'/'`;
    b=`echo $line | cut -f 2 -d'/'`;
		fn=`~/projects/PacBioSequencing/scripts/FindBasFiles.py $read --dir $2`
		outName=`echo $r | sed "s/:/_/"`
		outName="plots/$outName.$a.$b.pdf"
		echo "~/projects/PacBioSequencing/scripts/DotPlot.py --query bas:$fn:$b --target $3 --matches "dot:9" --savefig $outName --region $r"
		~/projects/PacBioSequencing/scripts/DotPlot.py --query bas:$fn:$b --target $3 --matches "dot:9" --savefig $outName --region $r
done < $1


