

    
    if (lchr not in allClusters):
        allClusters[lchr] = []
    clusters = allClusters[lchr]
    nClusters = len(clusters)
    ovpFound = False

    for i in range(nClusters):
        ovpRes = clusters[i].IntervalsOverlap(lchr, lstrand, lstart, lend, rchr, rstrand, rstart, rend)
        if (ovpRes == 1):
            clusters[i].Add(lstart, lend, 'leftBp')
            clusters[i].Add(rstart, rend, 'rightBp')
            ovpFound = True
            break
        if (ovpRes == 2):
            clusters[i].Add(lstart, lend, 'rightBp')
            clusters[i].Add(rstart, rend, 'leftBp')
            ovpFound = True
            break
    if (ovpFound == False):
        clusters.append(Cluster(lchr, lstrand, lstart, lend, rchr, rstrand, rstart, rend))
        if (nLines % 1000 == 0):
            sys.stderr.write( str(len(clusters)) + " clusters " + str(nLines) + " lines " + "\n")
            
        
    nLines += 1

for k in allClusters.keys():
    clusters = allClusters[k]
    for c in clusters:
        if (c.leftBp.Size() >= args.minSupport):
            leftBp = c.leftBp.Coordinates()
            rightBp = c.rightBp.Coordinates()
            if (leftBp[1] <= rightBp[0] or
                leftBp[0] >= rightBp[1]):
                outFile.write( str(c.leftBp.Size()) + "\t" + str(c.leftStrand) + "\t" + c.leftChr + "\t" + str(leftBp[0]) + "\t" + str(leftBp[1]) + "\t" + str(c.rightBp.Size()) + "\t" + str(c.rightStrand) + "\t" + c.rightChr  + "\t" +  str(rightBp[0]) + "\t" + str(rightBp[1]) + "\n")

if (outFile != sys.stdout):
    outFile.close()
