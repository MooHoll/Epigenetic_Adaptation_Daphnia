
##----------------------------------------------
# Eutrophic vs Recovery: DONE
##----------------------------------------------

file.list <- list("eutrophic_12_3_CpG.txt", "eutrophic_12_4_CpG.txt",	 "eutrophic_12_5_1_CpG.txt",  
                  "eutrophic_13_1_CpG.txt",	"eutrophic_13_2_CpG.txt",	"eutrophic_13_3_CpG.txt",	
                  "eutrophic_13_5_1_CpG.txt", "eutrophic_14_5_1_CpG.txt", "eutrophic_15_5_1_CpG.txt",
                  "recovery_0_1_CpG.txt","recovery_0_2_CpG.txt","recovery_0_4_CpG.txt",
                  "recovery_1_2_CpG.txt","recovery_2_1_CpG.txt","recovery_2_5_11_CpG.txt",
                  "recovery_2_5_9_CpG.txt", "recovery_3_5_15_CpG.txt","recovery_3_5_1_CpG.txt",
                  "recovery_3_5_2_CpG.txt", "recovery_3_6_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("12_3", "12_4", "12_5_1", "13_1", "13_2", "13_3",
                                      "13_5_1", "14_5_1", "15_5_1",
                                      "0_1", "0_2","0_4","1_2","2_1","2_5_11","2_5_9","3_5_15",
                                      "3_5_1","3_5_2","3_6"),
                     treatment = c(rep(0,9), rep(1,11)),
                     assembly="daphmag_2.4",
                     context="CpG")

##----------------------------------------------
# Eutrophic vs Pesticide: DONE
##----------------------------------------------

file.list <- list("eutrophic_12_3_CpG.txt", "eutrophic_12_4_CpG.txt",	 "eutrophic_12_5_1_CpG.txt",  
                  "eutrophic_13_1_CpG.txt",	"eutrophic_13_2_CpG.txt",	"eutrophic_13_3_CpG.txt",	
                  "eutrophic_13_5_1_CpG.txt", "eutrophic_14_5_1_CpG.txt", "eutrophic_15_5_1_CpG.txt",
                  "pesticide_6_2_CpG.txt","pesticide_6_3_CpG.txt","pesticide_7_3_CpG.txt",
                  "pesticide_7_5_CpG.txt","pesticide_7_5_4_CpG.txt","pesticide_8_5_3_CpG.txt",
                  "pesticide_9_5_1_CpG.txt", "pesticide_9_5_3_CpG.txt","pesticide_9_6_CpG.txt",
                  "pesticide_9_20_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("12_3", "12_4", "12_5_1", "13_1", "13_2", "13_3",
                                      "13_5_1", "14_5_1", "15_5_1",
                                      "6_2","6_3","7_3","7_5", "7_5_4","8_5_3","9_5_1",
                                      "9_5_3","9_6","9_20"),
                     treatment = c(rep(0,9), rep(1,10)),
                     assembly="daphmag_2.4",
                     context="CpG")

##----------------------------------------------
# Eutrophic vs Pristine: DONE
##----------------------------------------------

file.list <- list("eutrophic_12_3_CpG.txt", "eutrophic_12_4_CpG.txt",	 "eutrophic_12_5_1_CpG.txt",  
                  "eutrophic_13_1_CpG.txt",	"eutrophic_13_2_CpG.txt",	"eutrophic_13_3_CpG.txt",	
                  "eutrophic_13_5_1_CpG.txt", "eutrophic_14_5_1_CpG.txt", "eutrophic_15_5_1_CpG.txt",
                  "pristine_36_01_CpG.txt","pristine_36_02_CpG.txt","pristine_48_01_CpG.txt",
                  "pristine_48_02_CpG.txt","pristine_53_01_CpG.txt","pristine_54_01_CpG.txt",
                  "pristine_54_02_CpG.txt", "pristine_74_01_CpG.txt","pristine_77_01_CpG.txt",
                  "pristine_88_01_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("12_3", "12_4", "12_5_1", "13_1", "13_2", "13_3",
                                      "13_5_1", "14_5_1", "15_5_1",
                                      "36_01","36_02","48_01","48_02", "53_01","54_01","54_02",
                                      "74_01","77_01","88_01"),
                     treatment = c(rep(0,9), rep(1,10)),
                     assembly="daphmag_2.4",
                     context="CpG")

##----------------------------------------------
# Pristine vs Recovery: DONE
##----------------------------------------------

file.list <- list("pristine_36_01_CpG.txt","pristine_36_02_CpG.txt","pristine_48_01_CpG.txt",
                  "pristine_48_02_CpG.txt","pristine_53_01_CpG.txt","pristine_54_01_CpG.txt",
                  "pristine_54_02_CpG.txt", "pristine_74_01_CpG.txt","pristine_77_01_CpG.txt",
                  "pristine_88_01_CpG.txt",
                  "recovery_0_1_CpG.txt","recovery_0_2_CpG.txt","recovery_0_4_CpG.txt",
                  "recovery_1_2_CpG.txt","recovery_2_1_CpG.txt","recovery_2_5_11_CpG.txt",
                  "recovery_2_5_9_CpG.txt", "recovery_3_5_15_CpG.txt","recovery_3_5_1_CpG.txt",
                  "recovery_3_5_2_CpG.txt", "recovery_3_6_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("36_01","36_02","48_01","48_02", "53_01","54_01","54_02",
                                      "74_02","77_01","88_01",
                                      "0_1", "0_2","0_4","1_2","2_1","2_5_11","2_5_9","3_5_15",
                                      "3_5_1","3_5_2","3_6"),
                     treatment = c(rep(0,10), rep(1,11)),
                     assembly="daphmag_2.4",
                     context="CpG")

##----------------------------------------------
# Recovery vs Pesticide
##----------------------------------------------

file.list <- list("recovery_0_1_CpG.txt","recovery_0_2_CpG.txt","recovery_0_4_CpG.txt",
                  "recovery_1_2_CpG.txt","recovery_2_1_CpG.txt","recovery_2_5_11_CpG.txt",
                  "recovery_2_5_9_CpG.txt", "recovery_3_5_15_CpG.txt","recovery_3_5_1_CpG.txt",
                  "recovery_3_5_2_CpG.txt", "recovery_3_6_CpG.txt",
                  "pesticide_6_2_CpG.txt","pesticide_6_3_CpG.txt","pesticide_7_3_CpG.txt",
                  "pesticide_7_5_CpG.txt","pesticide_7_5_4_CpG.txt","pesticide_8_5_3_CpG.txt",
                  "pesticide_9_5_1_CpG.txt", "pesticide_9_5_3_CpG.txt","pesticide_9_6_CpG.txt",
                  "pesticide_9_20_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("0_1", "0_2","0_4","1_2","2_1","2_5_11","2_5_9","3_5_15",
                                      "3_5_1","3_5_2","3_6",
                                      "6_2","6_3","7_3","7_5", "7_5_4","8_5_3","9_5_1",
                                      "9_5_3","9_6","9_20"),
                     treatment = c(rep(0,11), rep(1,10)),
                     assembly="daphmag_2.4",
                     context="CpG")

##----------------------------------------------
# Pesticide vs Pristine: DONE
##----------------------------------------------

file.list <- list("pesticide_6_2_CpG.txt","pesticide_6_3_CpG.txt","pesticide_7_3_CpG.txt",
                  "pesticide_7_5_CpG.txt","pesticide_7_5_4_CpG.txt","pesticide_8_5_3_CpG.txt",
                  "pesticide_9_5_1_CpG.txt", "pesticide_9_5_3_CpG.txt","pesticide_9_6_CpG.txt",
                  "pesticide_9_20_CpG.txt",
                  "pristine_36_01_CpG.txt","pristine_36_02_CpG.txt","pristine_48_01_CpG.txt",
                  "pristine_48_02_CpG.txt","pristine_53_01_CpG.txt","pristine_54_01_CpG.txt",
                  "pristine_54_02_CpG.txt", "pristine_74_01_CpG.txt","pristine_77_01_CpG.txt",
                  "pristine_88_01_CpG.txt")

raw_data <- methRead(file.list,
                     sample.id = list("6_2","6_3","7_3","7_5", "7_5_4","8_5_3","9_5_1",
                                      "9_5_3","9_6","9_20",
                                      "36_01","36_02","48_01","48_02", "53_01","54_01","54_02",
                                      "74_02","77_01","88_01"),
                     treatment = c(rep(0,10), rep(1,10)),
                     assembly="daphmag_2.4",
                     context="CpG")