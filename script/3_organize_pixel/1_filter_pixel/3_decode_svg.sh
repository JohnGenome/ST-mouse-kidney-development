## Run in shell ##
# Extract coordinate
cd Filter_pixel_Annv5/origin
for si in E11_1 E11_2 E12_1 E12_2 E12_3 E12_4 E13_1 E13_2 E15_1 E15_2 E15_3 E15_4 E16_1 E18_1 E18_2 E18_3 E18_4
do
cat 1_spatial_map.${si}.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 1_spatial_map.${si}.coord
done
cd Filter_pixel_Annv5/outFinal
for si in E11_1 E11_2 E12_1 E12_2 E12_3 E12_4 E13_1 E13_2 E15_1 E15_2 E15_3 E15_4 E16_1 E18_1 E18_2 E18_3 E18_4
do
cat 3_spatial_map.${si}.final.svg | grep "path class" | sed 's/^.*M\([0-9]*\.\?[0-9]*,[0-9]*\.\?[0-9]*\).*$/\1/g' | sed 's/,/ /g' > 3_spatial_map.${si}.final.coord
done
