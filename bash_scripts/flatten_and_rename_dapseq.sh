# from https://unix.stackexchange.com/questions/45644/flatten-directory-but-preserve-directory-names-in-new-filename
#find . -type f -exec sh -c 'new=$(echo "{}" | tr "/" "-" | tr " " "_"); mv "{}" "$new"' \;

find . -depth -type f -exec rename 's@(?<!\.)/@_@g' -- {} \;

# remove empty directories
find . -type d -empty -delete

# clean filename 
for file in *.narrowPeak; do
    rename 's'
    mv "$file" "${file/_a_chr1-5_chr1-5_GEM_events/}"
    mv "$file" "${file/_b_chr1-5_chr1-5_GEM_events/}"
    mv "$file" "${file/_v3._chr1-5_chr1-5_GEM_events/}"
    mv "$file" "${file/_chr1-5_chr1-5_GEM_events/}"
done

for file in *.narrowPeak; do
    mv "$file" "${file/_tnt/}"
done