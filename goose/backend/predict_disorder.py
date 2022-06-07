#import stuff
import threading
from threading import Thread
import metapredict as meta

# Class to allow for return values with threading
class ThreadWithReturnValue(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs)
        self._return = None
    def run(self):
        #print(type(self._target))
        if self._target is not None:
            self._return = self._target(*self._args,
                                                **self._kwargs)
    def join(self, *args):
        Thread.join(self, *args)
        return self._return


# function for predicting disorder of last 20 aas of a sequence with a specific
# amino acid at the end
def predict(sequence, amino_acid, seq_range = []):

   """
   Function for predicting the disorder of a sequence either
   with a range or just the last 20 amino acids. Range is used
   for identifying regions in the sequence. Default is to just
   check the last 20 amino acids.

   Parameters
   -------------
   sequence : string
     the amino acid sequence as a string

   amino_acid : string
      The amino acid to add to the end of the sequence
      and to predict the disorder value for

   seq_range : list
      The coordinates of the substring to test
      default = empty list (range not used)

   Returns
   ---------

   Float
     Returns the predicted value for the disorder using the specified
     amino acid at the end of the sequence (or at the end of the sub sequence)
   """   

   # step 1. Figure out what the sequence is
   if seq_range == []:
      # cut off length to last 20 amino acids
      if len(sequence) > 20:
         used_sequence = sequence[len(sequence)-20:]
      else:
         used_sequence = sequence
   else:
      # get the subsequence
      sub_sequence = sequence[seq_range[0]: seq_range[1]]
      # cut off the length of the subsequence
      if len(sub_sequence) > 20:
         used_sequence = sub_sequence[len(sub_sequence)-20:]
      else:
         used_sequence = sub_sequence     
   
   # step 2. Add the amino acid to the sequence 
   final_sequence = used_sequence + amino_acid
   
   # step 3. Predict the disorder
   cur_prediction = meta.predict_disorder(final_sequence)

   # step 4. Grab the disorder prediction for the final vallue in the list which
   # corresponds to the value for the amino acid to be added to the sequence
   final_value = cur_prediction[-1]
   
   # step 5. return the final value
   return final_value



# function for calculating the disorder for a sequence using
# the predict function. Allows for simultaneous predictions.
def calculate_disorder_for_AA(sequence, amino_acid_list = [], seq_range=[]):
   
   """
   Function for predicting the disorder of a sequence either
   with a range or just the last 20 amino acids. Range is used
   for identifying regions in the sequence. Default is to just
   check the last 20 amino acids.

   Parameters
   -------------
   sequence : string
     the amino acid sequence as a string

   amino_acid_list : list
      The amino acids to add to the end of the sequence
      and to predict the disorder value for

   seq_range : list
      The coordinates of the substring to test
      default = empty list (range not used)

   Returns
   ---------

   Dict
     Returns the predicted value for the disorder using the specified
     amino acid at the end of the sequence (or at the end of the sub sequence)
     as a dictionary where the key:value pairs are the amino acid and the 
     predicted disorder for that amino acid after the input sequence.
   """     

   # if no amino acid list provided, use all
   if amino_acid_list==[]:
      amino_acid_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

   # create empty dict to hold results
   results_dict = {}

   # create thread(s)
   A_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'A', seq_range,))
   C_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'C', seq_range,))
   D_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'D', seq_range,))
   E_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'E', seq_range,))
   F_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'F', seq_range,))
   G_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'G', seq_range,))
   H_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'H', seq_range,))
   I_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'I', seq_range,))
   K_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'K', seq_range,))
   L_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'L', seq_range,))
   M_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'M', seq_range,))
   N_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'N', seq_range,))
   P_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'P', seq_range,))
   Q_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'Q', seq_range,))
   R_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'R', seq_range,))
   S_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'S', seq_range,))
   T_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'T', seq_range,))
   V_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'V', seq_range,))
   W_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'W', seq_range,))
   Y_thread = ThreadWithReturnValue(target=predict, args = (sequence, 'Y', seq_range,))

   # start necessary threads depending on which amino acids to test.
   # probably could have been done in a loop, but didn't want to mess with that.
   if 'A' in amino_acid_list: 
      A_thread.start()
   if 'C' in amino_acid_list: 
      C_thread.start()
   if 'D' in amino_acid_list: 
      D_thread.start()
   if 'E' in amino_acid_list: 
      E_thread.start()
   if 'F' in amino_acid_list: 
      F_thread.start()
   if 'G' in amino_acid_list: 
      G_thread.start()
   if 'H' in amino_acid_list: 
      H_thread.start()
   if 'I' in amino_acid_list: 
      I_thread.start()
   if 'K' in amino_acid_list: 
      K_thread.start()
   if 'L' in amino_acid_list: 
      L_thread.start()
   if 'M' in amino_acid_list: 
      M_thread.start()
   if 'N' in amino_acid_list: 
      N_thread.start()
   if 'P' in amino_acid_list: 
      P_thread.start()
   if 'Q' in amino_acid_list: 
      Q_thread.start()
   if 'R' in amino_acid_list: 
      R_thread.start()
   if 'S' in amino_acid_list: 
      S_thread.start()
   if 'T' in amino_acid_list: 
      T_thread.start()
   if 'V' in amino_acid_list: 
      V_thread.start()
   if 'W' in amino_acid_list: 
      W_thread.start()
   if 'Y' in amino_acid_list: 
      Y_thread.start()  

   # join necessary threads. Using the ThreadWithResults class, 
   # join() now returns a value, so we can use that to add the
   # amino acid / predicted disorder value pair to our results_dict.
   if 'A' in amino_acid_list: 
      results_dict['A']=A_thread.join()
   if 'C' in amino_acid_list: 
      results_dict['C']=C_thread.join()
   if 'D' in amino_acid_list: 
      results_dict['D']=D_thread.join()
   if 'E' in amino_acid_list: 
      results_dict['E']=E_thread.join()
   if 'F' in amino_acid_list: 
      results_dict['F']=F_thread.join()
   if 'G' in amino_acid_list: 
      results_dict['G']=G_thread.join()
   if 'H' in amino_acid_list: 
      results_dict['H']=H_thread.join()
   if 'I' in amino_acid_list: 
      results_dict['I']=I_thread.join()
   if 'K' in amino_acid_list: 
      results_dict['K']=K_thread.join()
   if 'L' in amino_acid_list: 
      results_dict['L']=L_thread.join()
   if 'M' in amino_acid_list: 
      results_dict['M']=M_thread.join()
   if 'N' in amino_acid_list: 
      results_dict['N']=N_thread.join()
   if 'P' in amino_acid_list: 
      results_dict['P']=P_thread.join()
   if 'Q' in amino_acid_list: 
      results_dict['Q']=Q_thread.join()
   if 'R' in amino_acid_list: 
      results_dict['R']=R_thread.join()
   if 'S' in amino_acid_list: 
      results_dict['S']=S_thread.join()
   if 'T' in amino_acid_list: 
      results_dict['T']=T_thread.join()
   if 'V' in amino_acid_list: 
      results_dict['V']=V_thread.join()
   if 'W' in amino_acid_list: 
      results_dict['W']=W_thread.join()
   if 'Y' in amino_acid_list: 
      results_dict['Y']=Y_thread.join()

   # return results dict
   return results_dict
      

