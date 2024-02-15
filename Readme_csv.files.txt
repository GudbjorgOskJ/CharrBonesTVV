This dataseries comes from a study on Arctic charr diversity.
Title: “Diversity in the internal functional feeding elements of sympatric morphs of Arctic charr (Salvelinus alpinus)”
Doi: https://doi.org/10.1101/2023.02.17.528955
Author: Guðbjörg Ósk Jónsdóttir (https://orcid.org/0009-0008-0502-5553), Laura-Marie von Elm, Finnur Ingimarsson (https://orcid.org/0000-0002-0815-7622), Samuel Tersigni, Sigurður Sveinn Snorrason, Arnar Pálsson (https://orcid.org/0000-0002-6525-8112), Sarah Elizabeth Steele (https://orcid.org/0000-0001-8404-5537). 
Status: Considered for publication in PLoS One

Readme about the csv files.
Below is information on all csv. files within this repository.
All files contain information on the same individuals however files do contain different subsets.
Master_sheet.csv, Articular_Angular_Ch.csv, Premaxilla_Ch.csv and Quadrate_Ch.csv have information on all 240 specimens.
External_Ch.csv has information on 208 specimens (missing: PC3831 to PC3902).
Dentary_Ch.csv has information on 239 specimens (missing: PC3846).
Maxilla_Ch.csv has information on 239 specimens (missing: PC4224).
Supramaxilla_Ch.csv has information on 232 specimens (missing: PC3928, PC3930, PC3933, PC4207, PC4209, PC4215, PC4224, PC4261).

Master_sheet.csv
Data sheet containing general information on all of the individuals used in this study.
Explanation of variables (columns). 
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
date_net_set = the date (D/M/Y) when nets were set into the lake.
date_net_up = the date (D/M/Y) when nets were taken out of the lake.
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
diphyll_inf = Diphyllobothrium infection, higher numbers more infected, lower numbers less infected.
age = the age (in years) of each individual.
logCsize_W = Log transformed centroid size of the whole body.
logCsize_H = Log transformed centroid size of the head.
logCsize_D = Log transformed centroid size of the dentary bone
logCsize_M = Log transformed centroid size of the maxilla bone
logCsize_A = Log transformed centroid size of the articular-angular bone
logCsize_Q = Log transformed centroid size of the quadrate bone
logCsize_P = Log transformed centroid size of the premaxilla bone
logCsize_S = Log transformed centroid size of the supramaxilla bone

Articular_Angular_Ch.csv
General information on each individual, only for individuals which had an undamaged Articular-angular bone used for this analysis.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Articular_Angular_LM.tps file.
new_image_name = Current names of images, as seen on Figshare
teeth_den = Average number of teeth for the dentary bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_max = Average number of teeth for the maxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pre = Average number of teeth for the premaxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pal = Average number of teeth for the palatine bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_vo = Teeth count for the vomer bone, pre individual.
teeth_glos = Teeth count for the glossohyal bone, pre individual.

Dentary_Ch.csv
General information on each individual, only for individuals which had an undamaged dentary bone used for this analysis.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Dentary_LM.tps file.
new_image_name = Current names of images, as seen on Figshare
teeth_den = Average number of teeth for the dentary bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_max = Average number of teeth for the maxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pre = Average number of teeth for the premaxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pal = Average number of teeth for the palatine bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_vo = Teeth count for the vomer bone, pre individual.
teeth_glos = Teeth count for the glossohyal bone, pre individual.

External_Ch.csv
General information on each individual, only for individuals which had external photographs taken of them.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Whole_Body_LM.tps and Head_LM.tps files.
new_image_name = Current names of images, as seen on Figshare

Maxilla_Ch.csv
General information on each individual, only for individuals which had an undamaged maxilla bone used for this analysis.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Maxilla_LM.tps file.
new_image_name = Current names of images, as seen on Figshare
teeth_den = Average number of teeth for the dentary bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_max = Average number of teeth for the maxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pre = Average number of teeth for the premaxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pal = Average number of teeth for the palatine bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_vo = Teeth count for the vomer bone, pre individual.
teeth_glos = Teeth count for the glossohyal bone, pre individual.
teeth_ang = Average angle of teeth for the maxilla bones, per individual. Left bone teeth angle and right bone teeth angle were combined and divided by two.

Premaxilla_Ch.csv
General information on each individual, only for individuals which had an undamaged premaxilla bone used for this analysis.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Premaxilla_LM.tps file.
new_image_name = Current names of images, as seen on Figshare
teeth_den = Average number of teeth for the dentary bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_max = Average number of teeth for the maxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pre = Average number of teeth for the premaxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pal = Average number of teeth for the palatine bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_vo = Teeth count for the vomer bone, pre individual.
teeth_glos = Teeth count for the glossohyal bone, pre individual.

Quadrate_Ch.csv
General information on each individual, only for individuals which had an undamaged quadrate bone used for this analysis.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Quadrate_LM.tps file.
new_image_name = Current names of images, as seen on Figshare
teeth_den = Average number of teeth for the dentary bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_max = Average number of teeth for the maxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pre = Average number of teeth for the premaxilla bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_pal = Average number of teeth for the palatine bones, per individual. Left bone teeth count and right bone teeth count were combined and divided by two.
teeth_vo = Teeth count for the vomer bone, pre individual.
teeth_glos = Teeth count for the glossohyal bone, pre individual.


Supramaxilla_Ch.csv
General information on each individual, only for individuals which had an undamaged supramaxilla bone used for this analysis.
num = the number of the individual in this particular sheet. The first individual listed is number one, and then so on.
id_me = id-tags used for this analysis. Each individual had their own id-tag.
id_Hafro = id-tags used by collaborates as some individuals were shared between this project and other projects.
id_marina = id-tags used by collaborates as some individuals were shared between this project and other projects. Note on some figures individuals are still tagged with this tagging scheme. 
species = factor variable. The species of each individual used in this analysis. Artic charr is labeled as B (Icelandic work for arctic charr is bleikja).
morph = factor variable. If individual was arctic charr they would be classified into one of the four charr eco-morph from lake Þingvallavatn. LB = large benthivorous charr, SB= small benthivorous charr, PL= planktivorous charr, PI = piscivorous charr.
length_cm = Fork length measurements (in centimeters, cm) for each individual.
weight_g = the weight of each individual, measured in grams, g.
sex = factor variable. The biological sex of each individual. F stands for female and M stands for male.
mat = Maturity index, on a scale from 0 to 7, where highest value were given to individuals with running gametes
age = the age (in years) of each individual.
location = name of the location, within lake Þingvallavatn, where nets (which caught this particular individual) were set.
old_image_name = Older name of images used during analysis. These older names are used in the Supramaxilla_LM.tps file.
new_image_name = Current names of images, as seen on Figshare
