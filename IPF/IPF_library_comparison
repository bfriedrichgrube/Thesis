> traml2 <- read.delim2('transitionlist_optimized.tsv', sep='\t', as.is=T)
> traml1 <- read.delim2('../spectrast2tsv.tsv', sep='\t', as.is=T)
> traml1_prottyp <- traml1[grep('^1/',traml1$UniprotID),]
> length(unique(traml1$UniprotID))
[1] 10608
> length(unique(traml1$UniprotID))
[1] 10608
> length(unique(traml1_prottyp$UniprotID))
[1] 9021
> length(unique(traml2$ProteinId))
[1] 8699
> length(unique(traml1_prottyp$PeptideSequence))
[1] 51988
> length(unique(traml2$PeptideSequence))
[1] 37002
> colnames(traml2)
 [1] "PrecursorMz"             "ProductMz"              
 [3] "PrecursorCharge"         "ProductCharge"          
 [5] "LibraryIntensity"        "NormalizedRetentionTime"
 [7] "PeptideSequence"         "ModifiedPeptideSequence"
 [9] "PeptideGroupLabel"       "LabelType"              
[11] "CompoundName"            "SumFormula"             
[13] "SMILES"                  "ProteinId"              
[15] "UniprotId"               "FragmentType"           
[17] "FragmentSeriesNumber"    "Annotation"             
[19] "CollisionEnergy"         "PrecursorIonMobility"   
[21] "TransitionGroupId"       "TransitionId"           
[23] "Decoy"                   "DetectingTransition"    
[25] "IdentifyingTransition"   "QuantifyingTransition"  
[27] "Peptidoforms"           
> length(unique(traml2$ModifiedPeptideSequence))
[1] 54449
> head(traml2$ModifiedPeptideSequence)
[1] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
[2] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
[3] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
[4] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
[5] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
[6] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
> length(grep('.(UniMod:1)',unique(traml2$ModifiedPeptideSequence))
+ )
[1] 563
> length(grep('(UniMod:7)',unique(traml2$ModifiedPeptideSequence))
+ )
[1] 14846
> head(unique(traml2$ModifiedPeptideSequence))
[1] ".(UniMod:1)AAAAKPN(UniMod:7)N(UniMod:7)LSLVVHGPGDLR"
[2] ".(UniMod:1)AAAAKPN(UniMod:7)NLSLVVHGPGDLR"          
[3] ".(UniMod:1)AAAAKPNN(UniMod:7)LSLVVHGPGDLR"          
[4] ".(UniMod:1)AAAAKPNNLSLVVHGPGDLR"                    
[5] ".(UniMod:1)AAAAVGAGHGAGGPGAASSSGGAR"                
[6] ".(UniMod:1)AAADGGGPGGASVGTEEDGGGVGHR"               
> unique(traml2$ModifiedPeptideSequence)[15:25]
 [1] ".(UniMod:1)AATASAGAGGIDGKPR"                                 
 [2] ".(UniMod:1)AC(UniMod:4)GLVASNLNLKPGEC(UniMod:4)LR"           
 [3] ".(UniMod:1)AC(UniMod:4)PLDQ(UniMod:7)AIGLLVAIFHK"            
 [4] ".(UniMod:1)AC(UniMod:4)PLDQAIGLLVAIFHK"                      
 [5] ".(UniMod:1)ADEDGEGIHPSAPHR"                                  
 [6] ".(UniMod:1)ADKPDM(UniMod:35)GEIASFDK"                        
 [7] ".(UniMod:1)ADKPDMGEIASFDK"                                   
 [8] ".(UniMod:1)ADVLDLHEAGGEDFAM(UniMod:35)DEDGDESIHK"            
 [9] ".(UniMod:1)ADVLDLHEAGGEDFAMDEDGDESIHK"                       
[10] ".(UniMod:1)AEAMDLGKDPN(UniMod:7)GPTHSSTLFVR"                 
[11] ".(UniMod:1)AEASSDPGAEEREELLGPTAQWSVEDEEEAVHEQC(UniMod:4)QHER"
> length(grep('(UniMod:)',unique(traml2$ModifiedPeptideSequence)))
[1] 27281
> length(grep('(UniMod:4)',unique(traml2$ModifiedPeptideSequence)))
[1] 12030
> length(grep('(UniMod:1)',unique(traml2$ModifiedPeptideSequence)))
[1] 563
> length(grep('.(UniMod:1)',unique(traml2$ModifiedPeptideSequence)))
[1] 563
> length(grep('K(UniMod:1)',unique(traml2$ModifiedPeptideSequence)))
[1] 0
> length(grep('K/(UniMod:1/)',unique(traml2$ModifiedPeptideSequence)))
[1] 0
> length(grep('(UniMod:35)',unique(traml2$ModifiedPeptideSequence)))
[1] 5427
> length(grep('K\(UniMod:1\)',unique(traml2$ModifiedPeptideSequence)))
Error: '\(' is an unrecognized escape in character string starting "'K\("
> length(grep('\(UniMod:1\)',unique(traml2$ModifiedPeptideSequence)))
Error: '\(' is an unrecognized escape in character string starting "'\("
> length(grep('\\(UniMod:1\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 563
> length(grep('.\\(UniMod:1\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 563
> length(grep('\\.\\(UniMod:1\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 220
> length(grep('K\\(UniMod:1\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 350
> length(grep('\\(UniMod:4\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 12030
> length(grep('\\(UniMod:7\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 14846
> length(grep('\\(UniMod:21\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 243
> length(grep('\\(UniMod:35\\)',unique(traml2$ModifiedPeptideSequence)))
[1] 5427
> length(unique(traml2$ModifiedPeptideSequence))
[1] 54449
> length(grep('\\(UniMod:',unique(traml2$ModifiedPeptideSequence)))
[1] 27281
> ls()
[1] "traml1"         "traml1_prottyp" "traml2"        
> a <- unique(traml1_prottyp$UniprotID)
> b <- unlist(strsplit(a, sep='/'))
Error in strsplit(a, sep = "/") : unused argument (sep = "/")
> b <- unlist(strsplit(a,'/'))
> d <- b[seq(2, length(b), by=2)]
> length(d)
[1] 9021
> length(intersect(d, traml2$ProteinId))
[1] 7624
> colnames(traml1_prottyp)
 [1] "PrecursorMz"           "ProductMz"             "Tr_recalibrated"      
 [4] "transition_name"       "CE"                    "LibraryIntensity"     
 [7] "transition_group_id"   "decoy"                 "PeptideSequence"      
[10] "ProteinName"           "Annotation"            "FullUniModPeptideName"
[13] "PrecursorCharge"       "PeptideGroupLabel"     "UniprotID"            
[16] "FragmentType"          "FragmentCharge"        "FragmentSeriesNumber" 
[19] "LabelType"            
> length(unique(traml1_prottyp$FullUniModPeptideName))
[1] 57671
> length(intersect(traml1_prottyp$FullUniModPeptideName, traml2$ModifiedPeptideSequence))
[1] 36846
> q()

