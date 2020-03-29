#Decompresses .gz files in a folder defined by ($1) and saves files to /unzipped/ folder

mkdir $1/unzipped

for file in $1/*.gz
do
zcat $file > $file.decomp
done

mv *.decomp $1/unzipped/
