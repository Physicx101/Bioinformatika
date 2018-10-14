#Menggunakan library biopython
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from Bio.HMM import MarkovModel
from Bio.HMM import Trainer
from Bio.HMM import Utilities


#profile HMM dari slide
    #AT--GA
    #A-C-CA
    #AC-AST

#Membuat array DNA dan state pada HMM
class arrayDNA(Alphabet.Alphabet):
    letters = ['A','C','G','S', 'T']

#I = Insert, D = Delete
class arrayState(Alphabet.Alphabet):
    letters = ['I', 'D', 'E','F', 'J','K','L','M','N','O', 'P']
    
#Membangung Hidden Markov Model
markovBuilder = MarkovModel.MarkovModelBuilder(arrayState(),arrayDNA())


#transisi dari semua state ke state lainnya
markovBuilder.allow_all_transition()




#probabilitas awal (berdasarkan profile, karena selalu mulai dari Match state maka M = 1)
markovBuilder.set_initial_probabilities({'M': 1, 'D': 0,'I':0})

#probabilitas untuk transisi dari satu state ke state lain
#state yang tidak diberi skor berarti tidak akan dilewati karena probabilitas 0
markovBuilder.set_transition_score('M', 'N', .67)
markovBuilder.set_transition_score('M', 'D', .33)
markovBuilder.set_transition_score('D', 'I', 1)
markovBuilder.set_transition_score('N', 'O', .5)
markovBuilder.set_transition_score('N', 'I', .5)
markovBuilder.set_transition_score('I', 'O', 1)
markovBuilder.set_transition_score('O', 'P', 1)

#emission probability dari state ke DNA
markovBuilder.set_emission_score('M', 'A', 1)
markovBuilder.set_emission_score('N', 'T', .5)
markovBuilder.set_emission_score('N', 'T', .5)
markovBuilder.set_emission_score('I', 'C', .5)
markovBuilder.set_emission_score('I', 'A', .5)
markovBuilder.set_emission_score('O', 'C', .33)
markovBuilder.set_emission_score('O', 'G', .33)
markovBuilder.set_emission_score('O', 'S', .33)
markovBuilder.set_emission_score('P', 'A', .67)
markovBuilder.set_emission_score('P', 'T', .33)


#Menginisialisasi Hidden Markov Model
markovModel = markovBuilder.get_markov_model()

#3 sequence yang akan dialign
seq1 = Seq('ATGA',arrayDNA())
seq2 = Seq('ACCA',arrayDNA())
seq3 = Seq('ACAST',arrayDNA())


#state untuk tiap sequence
seq1State = MutableSeq('MNOP',arrayState())
seq2State = MutableSeq('MDIOP',arrayState())
seq3State = MutableSeq('MNIOP',arrayState())


seq = [seq1, seq2, seq3]
states = [seq1State, seq2State, seq3State]

#training Hidden Markov Model dengan sequence di atas
trainer = Trainer.KnownStateTrainer(markovModel)
for i in range(len(seq)):
    trainingseq = Trainer.TrainingSequence(seq[i],states[i])
    trainedhmm = trainer.train([trainingseq])
    

    
#contoh query yang lain
testSeq = Seq('ATSA', arrayDNA())
testState = MutableSeq('MNOP', arrayState())

#mencari state terbaik untuk sequence dengan Viterbi Algorithm
predictedstates, prob = trainedhmm.viterbi(testseq, arrayState())

#mengeluarkan hasil probabilitas dari sequence, emission, statenya, dan predicted statenya 
print("Probabilitas Prediksi: %f" % prob)
Utilities.pretty_print_prediction(testSeq, testState, predicted_states)