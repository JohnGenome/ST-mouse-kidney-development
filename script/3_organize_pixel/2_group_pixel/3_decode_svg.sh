## Run in shell ##
# Extract coordinate
cd Plot/1_ai_edited_svg_E18
for si in E18_3
do
for li in L1 L2 L3 L4 L5
do
cat 2_spatial_map.${si}.${li}.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 2_spatial_map.${si}.${li}.coord
done
done
cat 3_spatial_map.E18_3.bound0.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 3_spatial_map.E18_3.bound0.coord
cat 3_spatial_map.E18_3.bound5.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 3_spatial_map.E18_3.bound5.coord
