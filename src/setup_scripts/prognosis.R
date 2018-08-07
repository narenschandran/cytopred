prognosis<- list()

prognosis["GSE10358"] <- list()
prognosis["GSE13159"] <- list() 
prognosis["GSE30285"] <- list()
prognosis["GSE61804"] <- list()

prognosis$GSE10358$good <- c("t(8;21)",
                             "t(8;21)(q22;q22)", 
                             "t(8;21)(q22;q22); -Y", 
                             "t(8;21)(q22;q22); del(9)(q12q22)",
                             "t(8;21)(q22;q22); del(9)(q21q31)",
                             "t(8;21)(q22;q22); del(9)(q22q34)",
                             "t(15;17)",
                             "t(15;17)(q22;q21)", 
                             "t(15;17)(q22;q21); + 8", 
                             "t(15;17)(q22;q21); trisomy 8", 
                             "t(15;17)(q22;q21); inv(12)(p12p13)", 
                             "t(15;17)(q22;q21); ider(17)(q10)t(15;17)(q22;q21)",
                             "inv(16)",
                             "inv(16)(p13q22)",
                             "inv(16)(p13q22); +Y", 
                             "inv(16)(p13q22); + 22", 
                             "t(16;16)(p13;q22); +22", 
                             "inv(16)(p13q22); +22", 
                             "inv(16)(p13q22); add(20)(q13.3)", 
                             "inv(16)(p13q22); t(2;4)(q34;q21)",
                             "del(7)(q22) trisomy 8  t(15;17)(q22;q21)",
                             "der(7)t(7;8)(p15;q13)t(8;21)(q22;q22)")

prognosis$GSE10358$intr <- c("Normal",
                             "trisomy 21",
                             "trisomy 8",
                             "del(9)(q21q33); del(20)(q12)", 
                             "trisomy 21; + 22",
                             "del(9q)",
                             "+8",
                             "t(9;11)(p22;q23)",
                             "trisomy 21  trisomy 3",
                             "trisomy 10; del(9)(q13q22)",
                             "add(7)(p22)",
                             "del(11)(p12)", 
                             "iso(11)(q10)",
                             "del(12)(p12.3p13.3)",
                             "add(9)(p22)",
                             "normal")

prognosis$GSE10358$poor <- c("del(5)(q22q33)",
                             "del(7)(q21q36)",
                             "t(9;22)(q34;q11.2)-7",
                             "t(11;19)(q23;p13.1)",
                             "del(5)(q22q33); t(10;11)(p13~p15;q22~q23) i(17)(q10)",
                             ">3",
                             "t(11;19)(q23;p13); del(20q)",
                             "+ 21(q22); del(7)(p11.2)",
                             "t(11;19)(q23;p13.1); inv(12)(p12p13)",
                             "t(11;19)(q23;p13); +X",
                             "t(3;3)(q21;q26); +6",
                             "trisomy 11", 
                             "trisomy 13; iso(7)(p10)",
                             "del(7)(p13p15)",
                             "t(8;10;21)(q22;q26;q22); -Y", 
                             "t(6;9)(p23;q34)",
                             "Poor Risk",
                             "complex [no favorable]",
                             "Complex",
                             "complex [incl")

prognosis$GSE10358$unknown <- c("N.D.",
                              "ND",
                              "ND",
                              "failed")


#prognosis$GSE10358$other <- c("-X",
#                              "trisomy Y",
#                              "-Xt(5;21)(q31;q22); t(1;6)(q21;p23)",
#                              "t(12;22)(p13;q12)",
#                              "-9  add(9)(p22)  +mar", 
#                              "8",
#                              "ND")

prognosis$GSE13159$good <- c("AML with inv(16)/t(16;16)",
                             "AML with t(15;17)",
                             "AML with t(8;21)")

# Be careful here. normal karyotype + other could technically include poor also

prognosis$GSE13159$intr_poor <- c("AML with t(11q23)/MLL",
                                  "AML with normal karyotype + other abnormalities")

prognosis$GSE13159$poor <-  c("AML complex aberrant karyotype")

prognosis$GSE13159$unknown <- c("ALL with hyperdiploid karyotype",
                              "ALL with t(1;19)",
                              "ALL with t(12;21)",
                              "c-ALL/Pre-B-ALL with t(9;22)",
                              "c-ALL/Pre-B-ALL without t(9;22)",
                              "CLL",
                              "CML",
                              "mature B-ALL with t(8;14)",
                              "MDS",
                              "Non-leukemia and healthy bone marrow",
                              "Pro-B-ALL with t(11q23)/MLL",
                              "T-ALL")


prognosis$GSE30285$good <- c("t(15;17)",
                             "inv(16)",
                             "t(8;21)")

prognosis$GSE30285$intr <- c("MLL-AF9",
                             "MLL-ENL")

prognosis$GSE61804$good <- c("AML with t(15;17), PML-RARA, FAB M3v/M3",
                             "AML with t(15;17), PML-RARA",
                             "AML with t(15;17), PML-RARA, FAB M3/M3v",
                             "AML with t(15;17), PML-RARA, FAB M3",
                             "AML with t(8;21), AML1-ETO",
                             "AML with t(8;21), AML1-ETO, therapy-related",
                             "AML with inv(16)/t(16;16), CBFB-MYH11, therapy-related",
                             "AML with inv(16)/t(16;16) kryptic rearrangement, CBFB-MYH11",
                             "AML with inv(16)/t(16;16), CBFB-MYH11")


prognosis$GSE61804$intr <- c("AML with normal karyotype",
                             "AML with normal karyotype, therapy-related",
                             "AML with normal karyotype, FLT3-LM(5)",
                             "AML with trisomy 13 sole",
                             "AML with trisomy 11 sole",
                             "AML with t(11q23)/MLL, t(9;11), MLL-AF9, therapy-related",
                             "AML with t(11q23)/MLL, t(9;11), MLL-AF9",
                             "AML with t(11q23)/MLL, t(9;11)",
                             "AML with t(11q23)/MLL, t(11;19), MLL-ENL, therapy-related",
                             "AML with t(11q23)/MLL, t(11;19)",
                             "AML with t(11q23)/MLL, t(11;19), MLL-ELL, therapy-related")


prognosis$GSE61804$poor <- c("AML with trisomy 8 sole",
                             "AML with complex aberrant karyotype, del(5q)",
                             "AML with complex aberrant karyotype, untypical",
                             "AML with complex aberrant karyotype",
                             "AML with complex aberrant karyotype, untypical, therapy-related",
                             "AML with complex aberrant karyotype, untypical, independent clones",
                             "AML with complex aberrant karyotype, del(5q), RAEB-2/AML",
                             "AML with inv(3)/t(3;3)",
                             "AML with del(5q) sole",
                             "AML with del(5q) + 1 additional aberration",
                             "AML with t(11q23)/MLL, t(6;11), MLL-AF6",
                             "AML with t(11q23)/MLL, t(10;11), MLL-AF10")

