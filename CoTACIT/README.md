1. The suffix of the raw sequencing file should be "_1.fq.gz" and "_2.fq.gz". 
2. Modify T7 and T5 barcodes in "01_T7_split.sh" and "03_T5_split.sh".
3. Copy all files here to the directory where the raw sequening files are located.
4. Run the command: bash oneForAll.spilt.sh
5. Run the same pipeline for mapping and deduplication as for TACIT.
