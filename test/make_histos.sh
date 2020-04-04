# ! the ntuple must have the following naming syntax: l1TNtuple-${sample}_${date}.root for this to work !
# Usage: ./make_histos.sh sample date
sample=$1
date=$2
echo "Compiling..."
g++ -g  -Wno-deprecated controlplot.C -o plot.exe -I$ROOTSYS/include -L$ROOTSYS/lib `root-config --cflags` `root-config --libs`
echo "Making histograms..."
for etarange in barrel endcap forward all; do
	for ptcut in 20 35 90 120 180; do
		./plot.exe l1TNtuple-${sample}_${date}.root plots_${sample}_${date}_l1_${ptcut}_${etarange}.root ${etarange} ${ptcut}
	done
done

exit 0
