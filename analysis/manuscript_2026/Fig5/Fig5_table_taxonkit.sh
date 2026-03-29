conda activate taxonkit

path_out=/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/output/250909_4474_samples/final_test/
db_path=/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/db/taxdmp

cd "$path_out"

cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d}" --lineage-field 2 --data-dir "$db_path" -F | taxonkit reformat -f "{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3,4 > ./taxids_lineage.txt
cat taxids.txt | taxonkit lineage --data-dir "$db_path" | taxonkit reformat -f "{d};{K};{p};{c};{o};{f};{g};{s}" --lineage-field 2 --data-dir "$db_path" -F | cut -f 1,3 > taxids_long.txt