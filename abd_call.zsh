#use java8 !!!!!!!!

#make fasta index
#smtools faidx ABD_identical.fa
#gatk-4.0.0.0/gatk CreateSequenceDictionary -R ABD_identical.fa

#CS
samtools merge temporal.bam ABD_identical/11-2YS_to_ABD_identical.bam ABD_identical/11-2YS_to_ABD_identical.bam
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I temporal.bam -O ABD_identical/CS_to_ABD_identical.addRG.bam --RGLB CS_library --RGPU machine --RGPL illumina --RGSM CS
samtools index ABD_identical/CS_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/CS_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/CS.vcf

#GC2
samtools merge temporal.bam ABD_identical/12-1PMI_to_ABD_identical.bam ABD_identical/12-1YS_to_ABD_identical.bam
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I temporal.bam -O ABD_identical/GC2_to_ABD_identical.addRG.bam --RGLB GC2_library --RGPU machine --RGPL illumina --RGSM GC2
samtools index ABD_identical/GC2_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/GC2_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/GC2.vcf


#GC2mut#1
samtools merge temporal.bam ABD_identical/13-2PMI_to_ABD_identical.bam ABD_identical/13-2YS_to_ABD_identical.bam
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I temporal.bam -O ABD_identical/mut1_to_ABD_identical.addRG.bam --RGLB mut1_library --RGPU machine --RGPL illumina --RGSM mut1
samtools index ABD_identical/mut1_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/mut1_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/mut1.vcf


#GC2 x CS
samtools merge temporal.bam ABD_identical/15-1YS_to_ABD_identical.bam ABD_identical/15-3PMI_to_ABD_identical.bam
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I temporal.bam -O ABD_identical/CSxGC2_to_ABD_identical.addRG.bam --RGLB CSxGC2_library --RGPU machine --RGPL illumina --RGSM CSxGC2
samtools index ABD_identical/CSxGC2_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/CSxGC2_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/CSxGC2.vcf


#GC2mut#2
samtools merge temporal.bam ABD_identical/25-3PMI_to_ABD_identical.bam ABD_identical/25-3YS_to_ABD_identical.bam
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I temporal.bam -O ABD_identical/mut2_to_ABD_identical.addRG.bam --RGLB mut2_library --RGPU machine --RGPL illumina --RGSM mut2
samtools index ABD_identical/mut2_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/mut2_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/mut2.vcf


#GC2mut#3
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I ABD_identical/26-3PMI_to_ABD_identical.bam -O ABD_identical/mut3_to_ABD_identical.addRG.bam --RGLB mut3_library --RGPU machine --RGPL illumina --RGSM mut3
samtools index ABD_identical/mut3_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/mut3_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/mut3.vcf

#GC2mut#3 x CS
gatk-4.0.0.0/gatk AddOrReplaceReadGroups -I ABD_identical/28-1YS_to_ABD_identical.bam -O ABD_identical/CSxmut3_to_ABD_identical.addRG.bam --RGLB CSxmut3_library --RGPU machine --RGPL illumina --RGSM CSxmut3
samtools index ABD_identical/CSxmut3_to_ABD_identical.addRG.bam
gatk-4.0.0.0/gatk HaplotypeCaller -I ABD_identical/CSxmut3_to_ABD_identical.addRG.bam -R ABD_identical.fa -O ABD_identical/CSxmut3.vcf
