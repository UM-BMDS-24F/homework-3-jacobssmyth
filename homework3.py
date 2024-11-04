
from Bio.Blast.Applications import NcbiblastnCommandline


blastn_cline = NcbiblastnCommandline(
    cmd="/Users/jake/Documents/SchoolStuff/Courses/649/ncbi-blast-2.10.1+/bin/blastp", 
    query="human.fa", 
    db="/Users/jake/Documents/SchoolStuff/Courses/649/homework3/mouse_db",    
    evalue=0.001, 
    out="blastttt.xml", 
    outfmt=5          
)

stdout, stderr = blastn_cline()

from Bio.Blast import NCBIXML

with open("blastttt.xml") as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    with open("blast_results.txt", "w") as out:
        out.write("Human_Sequence_ID\tMouse_ID\tAlignment\tE-value\tBit_Score\n")
        for blast_record in blast_records:
            human_sequence_id = blast_record.query
            if blast_record.alignments:
                best_alignment = blast_record.alignments[0]
                for hsp in best_alignment.hsps:
                    mouse_id = best_alignment.hit_id 
                    alignment = hsp.query[0:60] + "..." + hsp.match[0:60] + "..." + hsp.sbjct[0:60] 
                    e_value = hsp.expect  
                    bit_score = hsp.score  
                    out.write(f"{human_sequence_id}\t{mouse_id}\t{alignment}\t{e_value}\t{bit_score}\n")


#1
#In regards to which program was used, I chose blastp as it is typically used when dealing with protein sequences. 

#2
#When using blastp, the default substitution matrix is blosum62. According to the slides it's utilized for semi related sequences, a midrange. 

#3
#For the few parameters, 0.001 was used for detecting matrhces, so it's looking for more significant matches.
#Outfmt simply refers to xml. The others just refer to the files/db. 