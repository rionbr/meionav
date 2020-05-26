# Remove build
echo "-- Removing old build and output --"
rm -f *.class

# Build
echo "-- Build Main.java --"
javac -cp gephi-toolkit-0.9.2-all.jar PlotNetSingleMod.java

# Run
echo "-- Running all scripts --"
./run_enterocyte_dm.sh
./run_spermatocyte_dm.sh

./run_enterocyte_mm.sh
./run_spermatocyte_mm.sh

./run_enterocyte_hs.sh
./run_spermatocyte_hs.sh
