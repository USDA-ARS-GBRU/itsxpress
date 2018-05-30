import os

# This is the project Root
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))

# Define Taxa that Hmms were created for
# A = Alveolata
# B = Bryophyta
# C = Bacillariophyta
# D = Amoebozoa
# E = Euglenozoa
# F = Fungi
# G = Chlorophyta (green algae)
# H = Rhodophyta (red algae)
# I = Phaeophyceae (brown algae)
# L = Marchantiophyta (liverworts)
# M = Metazoa
# N = Microsporidia
# O = Oomycota
# P = Haptophyceae (prymnesiophytes)
# Q = Raphidophyceae
# R = Rhizaria
# S = Synurophyceae
# T = Tracheophyta (higher plants)
# U = Eustigmatophyceae
# X = Apusozoa
# Y = Parabasalia

taxa_choices = ["Alveolata", "Bryophyta", "Bacillariophyta","Amoebozoa", "Euglenozoa", "Fungi",
 "Chlorophyta","Rhodophyta","Phaeophyceae","Marchantiophyta","Metazoa","Microsporidia",
 "Oomycota","Haptophyceae", "Raphidophyceae"," Rhizaria","Synurophyceae",
 "Tracheophyta","Eustigmatophyceae","Apusozoa","Parabasalia"]

taxa_dict = {"Alveolata":"A.hmm","Bryophyta":"B.hmm", "Bacillariophyta":"C.hmm",
 "Amoebozoa":"D.hmm", "Euglenozoa":"E.hmm", "Fungi":"F.hmm","Chlorophyta":"G.hmm",
 "Rhodophyta":"H.hmm","Phaeophyceae":"I.hmm","Marchantiophyta":"L.hmm","Metazoa":"M.hmm",
 "Microsporidia": "N.hmm","Oomycota":"O.hmm","Haptophyceae":"P.hmm",
 "Raphidophyceae":"Q.hmm"," Rhizaria":"R.hmm","Synurophyceae":"S.hmm",
 "Tracheophyta":"T.hmm","Eustigmatophyceae":"U.hmm","Apusozoa":"X.hmm","Parabasalia":"Y.hmm"}