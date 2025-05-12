Generating sequence libraries in GOOSE
========================================
Apart from generating individual sequences, you can also generate libraries of synthetic sequences with predicted disorder that span user-defined ranges of sequence properties or fractions of amino acids. 
**NOTE**: This is likely the part of the codebase that is in the largest need of improvement. Apologies in advance!
I do plan on revamping this soon, but my current priorities are with updating core functionality. 

Generating sequence libraries by properties
--------------------------------------------

To generate sequence libraries by specifying properties, use the ``create.seq_property_library()`` function. An important things to note **GOOSE automatically gets rid of sequences not possible to make**. For example, a sequence with an NCPR value where the absolute value of the NCPR value is greater than the FCR will not be included because it's not possible to generate. In addition, GOOSE will check the hydropathy value and see if it is possible based on the FCR. If it's not possible, GOOSE will not make it. However **for incompatible hydropathy / charge values**, GOOSE will print out the sequences it was not able to generate and print them unless you set ``silent_failed_seqs`` to True. In this function, you can specify the length, FCR, NCPR, hydropathy, kappa, and disorder cutoff values. **Note** - GOOSE will do it's best to get the kappa value spot on, but it does allow for some error. It will also also adjust the FCR to match the NCPR if the two values are not compatible. This function will return a dictionary of sequences where each sequence is named after the property specified for generating the sequence.

 
For this function, you can specify the property as a single value, a list of two values, or a list with three values where GOOSE will make all possible values from the first number of the list to the second number of the list at an interval equal to the third value of the list.

**Example**

.. code-block:: python

    create.seq_property_library(40, FCR=0.1, NCPR=[0, 0.1], hydropathy=[2, 4, 0.5])

    {'>FCR_0.1_NCPR_0_hydropathy_2.0': 'QKPSQNKNHTPTGQGNSHPQDHPEQQQQQPPQQQSTQNTP', '>FCR_0.1_NCPR_0_hydropathy_2.5': 'NSNSTSENNKQNNGPHSPGTSQPPNFQPSAPQENSGNGKH', '>FCR_0.1_NCPR_0_hydropathy_3.0': 'SAPTQDPQSHYTQGGNEQTGGSPTGPPGWSHAKRSPSGQG', '>FCR_0.1_NCPR_0_hydropathy_3.5': 'NRSSGSCAPLNSAGGTTPGNKEVADPPPPGSTGSWGHQTH', '>FCR_0.1_NCPR_0_hydropathy_4': 'PSTHSSAGPSDTSASSSARSVPSSDSAVKSSCGSGASTTS', '>FCR_0.1_NCPR_0.1_hydropathy_2.0': 'QPPSPHQPLSHHSQQHNGNTKKQKSHQPNKNNSHPNNHNQ', '>FCR_0.1_NCPR_0.1_hydropathy_2.5': 'SGHSQGQNTHTKQGRQRGHGHVSPNQQHSSTPQHMQSPKT', '>FCR_0.1_NCPR_0.1_hydropathy_3.0': 'TSPSNHPQKPGPTPAGMQTGGTPGKTKHPHHPGSKLQQYT', '>FCR_0.1_NCPR_0.1_hydropathy_3.5': 'SAMLNASAGNPSGGQQRNSANLGPSRTTQKTSAQARSPTG', '>FCR_0.1_NCPR_0.1_hydropathy_4': 'PAPKPGAKVVSTSALQRVAKSSPPACSPGTHPGSSPTTSS'}

In the example above, the first value ``40`` in the function is the length. This value must be a single integer value and is **required**. The second value, ``FCR=0.1`` is specified as a single value, so all sequences generated will have that value. The third value ``NCPR=[0, 0.1]`` is specified as a list, so GOOSE made a set of sequences where NCPR was equal to the first value in the list and then a set where NCPR was equal to the second value in the list. Finally, ``hydropathy=[2, 4, 0.5]`` was set equal to a list with 3 values. This means GOOSE was told to generate sequences with a range of the first value of the list ``2``, to the second value in the list, ``4``, at an interval equal to the third value of the list ``0.5``. This resulted in the creation of sequences with values 2.0, 2.5, 3.0, 3.5, and 4.0.
**If the third value, which is the interval value, cannot be equally used between the range of sequences, GOOSE will just use the maximum value as the last value.** For example, hydropathy = [2, 4, 1.1] would result in hydropathy values of 2, 3.1, and 4.

**Additional Usage** - 
In addition, you can add a random name to the end of each sequence name by setting ``random_name=True``. This way, if you combine multiple libraries (if you want replicates of sequences with the varying properties you specify), you don't need to worry about overwriting anything due to a shared name. 

**Example**

.. code-block:: python

    create.seq_property_library(40, FCR=0.1, NCPR=[0, 0.1], hydropathy=[2, 4, 0.5], random_name=True)

    {'>FCR_0.1_NCPR_0_hydropathy_2.0_F0K7D2N5N6': 'TPNHTQPHKNHDNNNPSHNHTGNNPTNPQEHKGSQTNQPT', '>FCR_0.1_NCPR_0_hydropathy_2.5_G5R2A2L8F8': 'THQNPEDTHTTHPSMSRSNNPQLQNNGQRPAPPSSPHGHN', '>FCR_0.1_NCPR_0_hydropathy_3.0_E8V4W5C2N2': 'GPGSEHPHAPGDSSTGNNTSGPTKPSTGGALSQNRQPQYP', '>FCR_0.1_NCPR_0_hydropathy_3.5_C2C3D1F5P4': 'PPTQQPNGQSMSGGARHTTAAAAEGSARMAELSQNHSNGG', '>FCR_0.1_NCPR_0_hydropathy_4_I3S8I0G9W4': 'APTKGVAPETRSTSPAASSGAGGGGSSPASSMSPSSGDGS', '>FCR_0.1_NCPR_0.1_hydropathy_2.0_N8E1Y2R5D0': 'PPPSTGHQKQNYSQNHHNNPPQHQRWHRNGPPRPNSHQSG', '>FCR_0.1_NCPR_0.1_hydropathy_2.5_H3H7I2S0W7': 'NSGGSKRSSSPGPTPNQPQNGRNMPMPQNRQNHTNFQNTP', '>FCR_0.1_NCPR_0.1_hydropathy_3.0_I2K0M6E2T0': 'TRSHQQPQMHGMPSTSPNGCQTLNSPSSMRKGPPPQSGKN', '>FCR_0.1_NCPR_0.1_hydropathy_3.5_H3Q4E1L7W8': 'QTSSPQTMGRSQTTTGSASMQSSGMASTSRPPRFSSRSTG', '>FCR_0.1_NCPR_0.1_hydropathy_4_A4D7R8W6V7': 'RTSPTSVKPPTSACKTAAGSTPMTRSPSSSTLAVNGPPAP'}

Generating sequence libraries by fractions
-------------------------------------------

To generate sequence libraries by specifying fractions, use the ``create.seq_fractions_library()`` function. An important thing to note ***GOOSE automatically gets rid of sequences not possible to make***. This includes where the fraction of amino acids is greater than 1 (for obvious reasons) and if the fraction of any amino acid specified is greater than the max_aa_fractions limit. The max_aa_fractions limit can be manually overridden (see below). 

For this function, you can specify the amino as a single value, a list of two values, or a list with three values where GOOSE will make all possible values from the first number of the list to the second number of the list at an interval equal to the third value of the list. If only two values are given, GOOSE just uses those two values and does not assume any interval.

**Example**

.. code-block:: python

    create.seq_fractions_library(100, A=0.1, D=[0.1, 0.2], K=[0.1, 0.2, 0.05])

    {'>A_0.1_D_0.1_K_0.1': 'SKSSSGTKEHSEAEDNGEDGAATNHDHNDEHGRATGADDKKNHGHTKEHGAQHQSSQGNNNHDKSNSTRDAHNGARSDKRARNKEKQHQKGQAGENDHGE', '>A_0.1_D_0.1_K_0.15': 'EAAGGHQHGKRQSGKSQSADENKGRKKDESTNKDNTSRQSRETASQGKAKKNNGGPKKAGNQDDAQDESEGSQRSSQQAKDAKNGDDQTKTDEGHSTKAQ', '>A_0.1_D_0.1_K_0.2': 'TGGAGKDASAGDATKDRAKSDNKGKTKKERAAKTQNKSHNQAQEKRTGESSHKEKRKDGENQAKSHSKNHKQRKPADTQKTEDERHEEHGHEDKKDEDEQ', '>A_0.1_D_0.2_K_0.1': 'DSSNATTSNDQDDKSHDSNTQAHDQREVGNSKDSNSNASDDENKAGQENTSAEEDNPDDHEEKDDNDRDGHAKKKTSADKDGDNDREAKASHNGKNAEEG', '>A_0.1_D_0.2_K_0.15': 'DQGDEKAQTDASKSKNDTGAADKHGKAKQTGKEEENQDDGKTDDHSGPTDGQDNRGDKKSEGTDDKAKDQQDDTDEQATTTTKRAGHAADEDSTNRKKRS', '>A_0.1_D_0.2_K_0.2': 'SDDKDRRDKQAHSNKHADAKSNEASHRKKHAGKHGQDTGKKDDGQNKDSADKTHKTKDGDSEQKAHDTSEADQAKKDGDHNGEDGDEDGDKAGKKQGNKN'}

In the example above, the first value ``100`` in the function is the length. This value must be a single integer value and is **required**. The second value, ``A=0.1`` is specified as a single value, so all sequences generated will have that value. The third value ``D=[0.1, 0.2]`` is specified as a two item list, so GOOSE made a set of sequences where the fraction of **D** is either 0.1 or 0.2. Finally, ``K=[0.1, 0.2, 0.05]`` was set equal to a list with 3 values. This means GOOSE was told to generate sequences with a range of the first value of the list ``0.1``, to the second value in the list, ``0.2``, at an interval equal to the third value of the list ``0.05``. This resulted in the creation of sequences with values 0.1, 0.15, 0.2.
**If the third value, which is the interval value, cannot be equally used between the range of values from the lowest to the highest, GOOSE will just use start at the lowest value and increase until it can't any more and then will add the maximum value.** For example, ``K = [0.1, 0.2, 0.08]`` would result in K fraction values of , 0.1, 0.18, and 0.2.
In addition, you can add a random name to the end of each sequence name by setting ``random_name=True``. 

**Additional usage**

Some other things you can specify are:

``warn_user`` - used to determine whether to warn the user of any problems with generating the sequnces. By default is set to True, but you can set it to False. 
``robust_warnings`` - used to return a warning message FOR EVERY SEQUENCE THAT HAS A FRACTION VALUE NOT EQUAL TO THE INPUT VALUE. This can be annoying for large library generation and is by default set to False. Set to True to get more information on the sequences that have errors and what those errors are. 
``max_aa_fractions`` - used to override the max amino acid fractions. Input as a dict. Example below:

**Example**

.. code-block:: python

    create.seq_fractions_library(100, A=0.1, D=[0.1, 0.2], K=[0.1, 0.2, 0.05], max_aa_fractions= {'K': 0.16})

    {'>A_0.1_D_0.1_K_0.1': 'TQHHDEKNRRAEANDSPNGDEHAQDGKHSAEKQRTQHAENSDRDHSEGAKGNNQHGRKKQENERRAGGGKQTHKNTSQDHGDRNAKDDAQNGQQHHNKHA', '>A_0.1_D_0.1_K_0.15': 'QRKGNSANSGERADHTGDHDTQNAATTGKRKDNKEKKNDKHSARAQNTDKKAHTEKGSKHATQNAHNESQPGGDETNSKKHASTGQGKGNNDSRKGRDDN', '>A_0.1_D_0.2_K_0.1': 'AQQRQDQGGDAKDADDRTDARKDTETSPAKQEQAGRSDDKGPDDTDQKAESPTESNERDQQGQETGDDKQQKKGSEDAHSQDQDGKGPDQGKDAHQAGSR', '>A_0.1_D_0.2_K_0.15': 'EGPQATSTDDEDDHHKSKDESADEGAKSGKRTEENDRAATDTAHHATKDHDDHHKTDGPEKTKDETKKADKEGHHKKDTAEKQEDNANSSTDDTQPSKDD'}

In the above example, we manually overrode the max fraction for K and set it to 0.16. This eliminated sequences where the K fraction was 0.2 like in the example above where the max_aa_Fractions were left as default.
