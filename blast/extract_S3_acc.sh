blastdbcmd -db "/shared/db/genbank/256/blast/nr" -entry all -outfmt "%a %t" | grep -E "(rps3 |ribosomal protein S3 )" | awk '{print $1}' >> acc_list.txt

blastdb_aliastool -db "/shared/db/genbank/256/blast/nr" -seqidlist rps3_acc.txt -dbtype prot -out rps3