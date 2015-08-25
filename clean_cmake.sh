# This script is just intended to make the development
# of ParMooN's CMake system easier. It is not supposed
# to be used unauthorized.
# CB 2015/08/20

echo "Step back. I'm cleaning ParMooN CMake."
rm CMakeCache.txt
rm -rf CMakeFiles
rm -rf Examples/CMakeFiles
rm -rf src/AMG/CMakeFiles
rm -rf src/FE/CMakeFiles
rm -rf src/General/CMakeFiles
rm -rf src/Geometry/CMakeFiles
rm -rf src/Parallel/CMakeFiles
rm -rf src/QuadFormulas/CMakeFiles
rm -rf src/Refinement/CMakeFiles
rm -rf src/System/CMakeFiles
echo "I'm finished."