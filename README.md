![GOOSE_LOGO_FULL](https://github.com/ryanemenecker/goose/blob/main/images/goose_logo_3.png) 
# GOOSE : Generate disOrdered prOtiens Specifying paramEters

## What is GOOSE?
GOOSE is a python package developed to make generating IDRs with specific properties or IDR variants easy. This was the original goal, but some of the functionality has expanded into other areas such as generating sequences that have predicted secondary structure. Basically, we want to make it easy for you to design IDRs or IDR variants for your research. GOOSE currently only works from Python, but a CLI is in the works. 

## WARNING - The current release is a 0.1 release. There are possible remaining bugs that we are working on, so it is good practice to double check that your generated sequence was what you intended. You can easily do this using the analyze functionality!

In addition, some upcoming features are still being finalized such as library creation, a CLI, and prediction of NES (the current implementation is sub-optimal as more data is needed). We are always working on improving our protein feature predictors, so you can expect that to improve in the near future, too.

We are also very excited to hear from you about any feature requests! Please email me at remenecker@wustl.edu with the subject line **GOOSE feature request** and I'll be happy to see if I can make it happen. My goal is to make this *as useful as possible to the IDR community*.

## Why was GOOSE made?

It is rather difficult to know which properties of any given intrinsically disordered region (IDR) are important for its function. Is it the hydropathy? Maybe the fraction of charged residues? Who knows? Further, generating IDR variants can be a tricky task. Sure, you can swap all of one residue for another or make large truncations in your IDR, but there are limitations to those and similar approaches. This is why we made GOOSE. GOOSE lets you generate IDRs or IDR variants based on an IDR that you are interested in to figure out which properties of your IDR of interest are important for its function. 

## What all can GOOSE do?

There are four main functionalities currently in GOOSE. 

**1.** GOOSE can let you generatae IDRs where you can specify the paramters *length*, *average hydropathy*, *fraction of charged residues (FCR)*, *net charge per residue (NCPR)*, and *kappa*, which is a parameter that defines charge separation in a sequence. These parameters are all very important for the biophysical behavior of any given IDR! 

**2.** GOOSE can let you generate IDR variants. The different variants I currently have implemented should allow you to address many questions you might be interested in by letting you choose which parts of your IDR to hold constant. Importantly, most of the variant generators **keep the overall properties of your idea constant**, unles you specifically decide to change it. This means that you can generate numerous sequences with the same overall properties but differing in slightly different ways to really narrow down on *what makes your IDR of interest work*?  More on that later though.

**3.** GOOSE can let you generate sequences **predicted** to form secondary structurers including alpha helices, beta strands, and beta sheets. Wait! This is for IDRs?! Fair point! However, I thought that because IDRs can interact with folded domains, this might be a feature that allows researchers to carefully investigate how their IDR of interest is impacted by the presence of, for example, a helix-turn-helix or a large beta sheet. I'm also working on making it so you can control some parameters of these to further refine how you can investigate the interactions between folded domain properties and your favorite IDR, but that's proven to be a bit challenging and will be coming down the line.

**4.** GOOSE comes with a bunch of sequence analysis tools. This is for a few reasons. First, sometimes it's nice to just analyze an IDR you're interested in before making variants. It's good to know whether your IDR may have a transcripitional activation domain or a phosphosite. Additionally, it allows you to analyze the IDRs you generate using GOOSE to make sure you don't accidentally sent it to the nucleus or something of the sort due to it having an inadvertantly added NLS.

## How can I use GOOSE?

**GOOSE is currently a Python API only.** I will be bringing out the command-line interface soon! Just focusing on feature implementation and bug fixes / extensive testing for the time being. 

## Installing GOOSE - GOOSE takes flight!

Right now to install GOOSE, you need to install it through Github 

To clone the GitHub repository and gain the ability to modify a local copy of the code, run

	$ git clone https://github.com/ryanemenecker/goose.git
	$ cd goose
	$ pip install .

I will try to get GOOSE up on PyPi ASAP so it is easier to install!

**Important note**

*AFTER you install GOOSE*, you will also need to install SPARROW (another package we ware working on that is not yet on Github). I promise we will make this easier in the near future, but for now you will also need to do the following:

	$ pip install git+https://github.com/holehouse-lab/sparrow.git

which will install SPARROW, a package necessary for my being able to generate sequences with specified kappa values. **Important note**: if your attempted install of SPARROW fails, it may be because you do not have numpy or cython installed. I made them both required for installation of GOOSE, so if you install GOOSE first, you should be ok. If you still have problems, try running from terminal:

	$ pip install cython
	$ pip install numpy

Sorry this is a pain, we promise we will make it better shortly!!!

# Time to leave the nest - how to use GOOSE from Python

## Sequence Generation

GOOSE allows you to easily generate sequences. The way this works for sequence generation is that GOOSE will first make a sequence that is *likely* to be disordered that has your desired sequence parameters. GOOSE will then verify the disorder of the generated sequence using metapredict V2. Once verified to be disordered, the sequence is returned back to you. The default cutoff disorder is 0.6, which is fairly strict. However, you can alter this if you'd like. For example, if you're having a hard time making a sequence, it may be worth reducing the threshold to the normal metapredict V2 threshold of 0.5.

## Some important to know limitations

**A quick note on discrepencies between input parameters and the actual values of generated sequences**. GOOSE by default allows a small amount of error for hydropathy (it is allowed to be of by 0.05, which is honestly negligible) and for kappa (it is allowed to be off by 0.02). This is basically a balance of accuracy and speed. The more constraints on a sequence, the more precise GOOSE hase to be, the longer it takes. If you go into goose/backend/parameters, you can modify these parameters globally if you install GOOSE locally. Finally, if you choose an NCPR / FCR combination that is mathematically impossible, GOOSE will get as close as it can but despite how feirce GEESE are, even GOOSE has to let math win that round.


### An important note on the speed of sequence generation

The protein disorder field is not a slow moving one, and we are not here to slow you and your research down. Therefore, it was important for us to make GOOSE as fast as possible. However, because GOOSE incorporates stochasticity into sequence generation (see below for more on that), GOOSE still has to do some work when designing your disordered sequence. Importantly *the more parameters you specify, or the more constraints you put on sequence design, the more time it will take GOOSE to generate your sequence*. 

 
### Stochasticity of sequences generated using GOOSE

Making a sequence generator to generate IDRs is not a challenging problem provided you simply choose the most disordered residue possible at each step. However, that would result in the same sequence being generated over and over again, which massively limits the utility of the sequence generator. We designed GOOSE to have a bit of stochasticity in that it will *generally* not make the same sequence twice. For example, when we made 100,000 sequences that were 100 amino acids in length using GOOSE, every single sequence was unique. However, we want to emphasize that *the more parameters you specify when making sequences, the less available sequence space there is for GOOSE to explore. This will increase the chance that GOOSE will make the same sequence twice*. This is generally not a problem, but it is something that we wanted users to be aware of.

*An important exception to this rule is the ``create.minimal_var()`` function because it attempts to change as little of the input sequence as possible and therefore will frequently return the same sequenece multiple times.*


### An important note on failed sequnce generation

Sometimes GOOSE can't make your sequence. We did our best to make sure that you will not run across a situation where this occurs, but with the enormous amount of possible different combinations of parameters and sequence properties, it is possible that GOOSE will be unable to make something you would like. This can often be easily overcome by *slightly* adjusting your parameters or reducing the cutoff value (the cutoff value is the value required for something to be considered disordered).


### Incompatible FCR, NCPR, and hydropathy values**

GOOSE uses a rescaled Kyte Doolittle hydropathy scale for calculating mean hydropathy. This scale goes from 0 to 9 where higher values are more hydrophobic.  Importantly, each amino acid has a specific value associated with it. For example, isoleucine (I) has a value of 9 and arginine (R) has a value of 0. The charged residues have low hydropathy values (R = 0, K = 0.6, D = 1, E = 1). Therefore, if you have a sequence with too many charged residues, you limit how high the mean hydropathy can possibly go for a sequence. If you specify a high FCR (or a low or high NCPR) and a high hydropathy, that sequence may be mathematically impossible to make. GOOSE will return an error if you do this. In addition, there is a limited sequence space of charge / hydropathy combinations that are predicted to be disordered. I have empirically determined the limits of FCR values associated with hydropathy values for GOOSE, and GOOSE will not let you cross that limit when designing sequences. 


### Limits on mean hydropathy**

The more hydrophobic a sequence, the less likely it is to be disordered *generally speaking*. The higher the mean hydropathy value, the harder time GOSE will have in making it.

**Importantly, all of these limits are also in place for when you are generating sequence variants!**


# Using GOOSE from within Python

To use GOOSE from within Python, first import create from goose

    from goose import create

Once GOOSE has been imported, you can start making seqeunces and sequence variants!

## Generating sequences with specified sequence parameters

The ``create.sequence()`` function lets you create sequences predicted to be disordered with various specified properties. With this function you must specify length (first arguement) and can also specify: 1. hydropathy, 2. fraction of charged residues (FCR), 3. net charge per residue (NCPR), and 4. kappa (charge asymmetry paramter where higher values mean more chager asymmetry). *5.* **sigma (can only specify sigma alone, can't combine with other properties).**

Here is some more info on the various arguemnts - 

``FCR`` - The fraction of charged residues in the sequence.
    You can also input ``fcr`` or ``fraction`` when specifying this argument
        Values must be between 0 and 1.

``NCPR`` - The net charge per residue of the sequence.
    You can also input ``ncpr`` or ``net_charge`` when specifying this argument
        Values must be between -1 and 1

``hydropathy`` - The mean hydropathy of the sequence.
    You can also input ``mean_hydro``, ``Hydro``, or ``hydro`` when specifying this 	argument. Values must be between 0.0 and 6.1.

``kappa`` - The kappa value of the sequence. 1 is maximum charge asymmetry, 0 is 		minimal asymmetry. The value must be between 0 and 1. The function can have a 			hard time hitting specific values if there are few charged residues, so you 			might have to adjust this slightly when making sequencees

In addition, it's important to note that ``sigma`` values must be between 0 and 1.


Finally, I added in a 'spell check' function, so if you accidently spell something wrong like *hydropath* instead of *hydropathy*, it'll still work for ya! 

**Importantly, you do not need to specify all of these sequence parameters simultaneously.** For example, if you specify FCR and hydropathy, GOOSE will return sequences that have varying NCPR values while making sure that the specified hydropathy and FCR values are what you input. In this way, you can generate many sequences that have fixed parameters that you want to stay fixed while other parameters can vary randomly. Alternatively, if you need to specify values for all parameters, you can do that too!

**Examples**

Just specifying sequence length -

    create.sequence(40)
    'GDHNKAGQPPRKCSDQGGAGAPNPDCDPDTAPMDGDRMTN'

**Note** the length must between 10 and 10,000!

**Specifying additional properties -**

**Hydropathy**

    create.sequence(100, hydro = 3)
    'MTSYGRDGSPETGEGSTGTNSSSSRSMMGSTHNWQQYNGGTTSGTSSTGDSHRTHGDHSAGETTSGGDSEGTDETSTTTNGRGSSSGHDGSTGQDTNTRR'


Hydropathy values can be between 0.0 and 6.1. **Note**: the higher the hydropathy *over 5*, the longer it will take GOOSE to generate the sequence. Sequences that are very hydrophobic and disordered can be tricky to make. **Note**: whenever specifying hydropathy values, GOOSE will return a sequence within 0.05 of the specified value! This amount of error helps keep GOOSE fast (and a difference of 0.05 is typically negligible).


**fraction of charged residues (FCR)**

    create.sequence(40, FCR = 0.3)
    'GDRPSEHGQGPRKEDGMDQDDVSTEGHEWSNNPCNQSNNP'

FCR values can be between 0 and 1

**net charge per residue (NCPR)**

    create.sequence(40, NCPR = -0.2)
    'MQKNDRAPDHKDREKDGPIKERPEECPDDEQSDDEECPSH'

NCPR values can be between -1 and 1.

**Sigma**

    create.sequence(40, sigma = 0.3)
    'EKDKMEETHDDEGMQQDNNTETDEQPDNYESNDDEHATEG'

NCPR values can be between 0 and 1.

 
**Specifying multiple properties**

GOOSE lets you combine different properties simultaneously. Possible combinations are: FCR, NCPR, kappa & hydropathy can all be specified at once if you want! Importantly, any value you do not specify will just be random.

**Note**: Specifying high hyrdropathy values that have high FCR are challenging to make and may be slow! In addition, simply because of the hydropathy value of charged resides (K=0.6, R = 0, D and E = 1), it is impossible to make sequences with too high of fraction of charged residues and high hydropathy values. 

**Examples**

**FCR & NCPR**

    create.sequence(100, FCR = 0.3, NCPR = -0.1)
    'TSNQDKEMPQQHSPRCQPGEKVSDPPRSSDNSTNGGARPQQDWRPPEHMNPNRYEPNTMHQNREGRESAGGKDWPNPTIDQNQDPHEDTDNQEEESDHPC'

You cannot have values for NCPR where the absolute value of NCPR is greater than the specified FCR value. For example, NCPR = 0.4, FCR = 0.2 will not work (you can't get a net charge of 0.4 with only 0.2 fraction of charged residues!)

**Important note on combining FCR and NCPR!** Whenever NCPR and FCR are combined, if the combinations of the length, NCPR, and FCR is not mathematically possible, GOOSE will get as close as it can. In addition, GOOSE will prioritize NCPR over FCR and may have to slightly change the FCR to keep the NCPR accurate. This was necessary in that one of the 2 had to be prioritized, and it made the most sense to prioritize net charge.


**FCR & Hydropathy**

    create.sequence(100, FCR = 0.3, hydro = 3.2)
    'KVDSGTTSCSGERESDSGDLKSSKEGSSGSGSSSKSSKSKEATGSSTDTTAAAGGKGGGGGGDGGKGDGRGKGGGGGGEGRDGGGGGGEGGRGGGGRKRD'

When specifying hydropathy with FCR or NCPR, the max possible hydropathy value is 5.2. In addition, after extensive testing I found that sequences with high hydropathy values and high FCR values will never be predicted to be disordered by metapredict. Therefore, I restricted the ability to input these parameter combinations (no sense being stuck waiting for a sequence to be generated that will never actually end up generated). In general, the maximum possible FCR value will equal (hydro x -0.2289)+1.2756.

**NCPR & Hydropathy**

    create.sequence(100, NCPR = -0.3, hydro = 2.4)
    'REARGDAKGERDRGGDAKDKGAESGKDDDGEEEGAGEEEGEEGDDEAEADRADKERAERDKGDRDRAEGRAEKGAAAAEGADEGADEADEEEDDDADDEE

When specifying hydropathy with FCR or NCPR, the max possible hydropathy value is 5.2.


**NCPR, FCR, & Hydropathy**

    create.sequence(100, hydro = 2.65, NCPR = -0.3, FCR = 0.4)
    'NETPARPETHRDTASTSEGDETSEPEGTWSSNEADTDDDAETEHSPMSEDGERCESSKDAPPMRDEEGDDEDVEDTPDVSSSPDYEPGGHYSESNNDWPD'

This function has the same limitations as FCR & hydropathy or NCPR & hydropathy. It's important to note that **the more properties you specify, the longer it will take GOOSE to make the sequence**. In addition, longer sequences tend to take greater amounts of time to generate.


**NCPR, FCR, Hydropathy, and kappa**

    create.sequence(100, hydro = 2.65, NCPR = 0.0, FCR = 0.4, kappa=0.2)
    'NETPARPETHRDTASTSEGDETSEPEGTWSSNEADTDDDAETEHSPMSEDGERCESSKDAPPMRDEEGDDEDVEDTPDVSSSPDYEPGGHYSESNNDWPD'


### Generating sequences by specifying fractions of amino acids in Python

The ``create.seq_fractions()`` function lets you create sequences predicted to be disordered with specified fractions of various amino acids. With this function, you can specify multiple amino acids simultaneously, and each fraction should be specified using a decimal value (for example, if you want one tenth of the amino acids to be alanine use A=0.1).

For each amino acid, possible maximum values are as follows - 

"A" - 0 : 0.9, 
"R" - 0 : 1.0, 
"N" - 0 : 1.0, 
"D" - 0 : 1.0, 
"C" - 0 : 0.16, 
"Q" - 0 : 0.72, 
"E" - 0 : 1.0, 
"G" - 0 : 1.0, 
"H" - 0 : 1.0, 
"I" - 0 : 0.2, 
"L" - 0 : 0.26, 
"K" - 0 : 1.0, 
"M" - 0 : 0.26, 
"F" - 0 : 0.18, 
"P" - 0 : 0.94, 
"S" - 0 : 0.88, 
"T" - 0 : 0.76, 
"W" - 0 : 0.22, 
"Y" - 0 : 0.22, 
"V" - 0 : 0.3

**Examples**

**Just specifying a single amino acid -**

    create.seq_by_fractions(100, Q=0.3)
    'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQSTKMHDQHDSVNEDQEQNQNPWGHQPHMKGESNSSAREAQSEDQQNQAQNQQQNHDSTQQQDGQMDQ'

**Note**: Some fractions of amino acids are simply not possible. For example, GOOSE cannot make a sequence with W (tryptophan) greater than 0.2!


**Specifying multiple amino acids -**

    create.seq_by_fractions(100, Q=0.3, S=0.3, E=0.1)
    'QEQQSQKASQSQVESQDSSESSAPGSSQMHQQQSQSQEGMEQHQSSVGNSSSYPQSEQSEQQRQQSSQDQQQQSSSQTSEENSQSRQHDMSDTEMSGSQR'

**Note** - 
Some combinations of amino acids are simply not possible to make disorded. For example, while you can have up to 0.2 for W, if you have 0.18 for W and 0.16 for Y (both below their maximums), the total number of aromatic residues makes it unlikely that the sequence will be predicted to be disordered.


**Specifying none of a specific amino acids -**

    create.seq_by_fractions(50, A=0)
    'NKERPTGSWDEPPFDEGSSGMTNEDMGNKPYPTTDMQPEKWPQNDQQGST'

 

## Creating Sequence Variants in Python

Apart from simply generating sequences, GOOSE can help you make different types of sequence variants. In contrast to when you generate sequences, the primary input into the sequence variant functions is your sequence of interest. This makes things a bit more streamlined in that you don't have to figure out the properties of your sequence variant beforehand, GOOSE will take care of that for you. Of course, if you'd like to know the properties of your sequence of interest before making variants, GOOSE can help you with that too. 

### An important note on disorder cutoffs when creating sequence variants

One problem we encountered when finding ways for people to easily make sequence variants occurs when something with a low predicted disorder score is used as the input sequence. In this situation, if GOOSE has a hard line cutoff disorder value that it must get to in order to return the sequence variant to you, it will almost certainly be unable to do so. To bypass this issue, GOOSE will first examine the disorder of your input sequence that is being used for sequence variant generation. Then, GOOSE will use the disorder values across your input sequence as the threshold values for disorder predictions. In this way, GOOSE will not necessarily create something more disordered than your input sequence. Nonetheless, all returned variants will have at a minimum the same amount of predicted disorder as the iput sequence *as predicted by metapredict*.


#### Types of sequence variants

``new_seq_constant_class_var()`` - A function to generate a variant where the sequence composition is new but the numbers of each residue from each class is the same. The overall properties of the generated sequence will also be constant.

``constant_class_var()`` - function to generate a variant with the same properties as the input variant as well as the same order of amino acids as far as class and the same number in each class. It will try to change the sequence as much as possible within these restraints.

``new_var()`` - Function to generate a variant that is completely different in sequence to the input but has all the same overall parameters. Does not account for specific classes of residues.

``constant_residue_var()`` - function that will generate a new sequence variant where specific residues are held constant. The variant will have the same aggregate properties as the original sequence.

``shuffle_var()`` - Variant that will shuffle specific regions of an IDR. Multiple regions can be specified simultaneously.

``constant_class_hydro_var()`` - Function to take in a sequence and makes a variant that adjusts the hydropathy while keeping the position and number of amino acids the same by class of amino acid. Hydropathy is simply adjusted by changing what the residue at that position is by swapping that residue with  others that are in the same class but have differing hydropathies.

``kappa_var()`` - Variant where you can alter the charge asymmetry by changing the kappa value. Requires the presence of positively charged and negatively charged residues in the original sequence. Higher kappa values increase charge asymmetry, lower kappa values reduce charge asymmetry. Values can be between 0 and 1. 

``minimal_var()`` - Function for generating a variant that will change any parameters including hydropathy, fcr, ncpr, and scd (sequence charge decoration, which is another measurementof sequence charge asymmetry) to those that you want to change while minimizing the number of residues changed in the returned variant from the original sequence. The more you change a given parameter, the more the returned variant sequence will differ from the original.

For all class variants, the classes are categorized as followed:

aromatic : 'F', 'W', 'Y' 
polar : 'Q', 'N', 'S', 'T' 
positive : 'K', 'R' 
negative : 'D', 'E' 
hydrophobic' : 'I', 'V', 'L', 'A', 'M'
Special Cases : 'C', 'P', 'G', and 'H'
The 'Special Cases' residues are, for any function that accounts for the class of a residue, not interchangable with any other residues. 

### The new_seq_constant_class_var()

The ``new_seq_constant_class_var()`` keeps the same number of each class of residues but will return a sequence with a new overall sequence composition The overall properties of the generated sequence will also be constant.

**Example**

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.new_seq_constant_class_var(test)
    NENAPEDQENNNRNPNGNQNQANQNGPFDSAGT


### The constant_class_var()

The ``constant_class_var()`` generates a variant that keeps the order of amino acids the same by class but changes them as much as possible **while keeping the overall properties of the returned variant the same as the original sequence**. 

**Example**

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_class_var(test)
    NDNNGAENNDQQPRNEYPGNNNPTNNADGSNAS


### The new_var()

The ``new_var()`` generates a variant that is completely different in sequence to the input. It does not account for classes of residues. The returned sequence will simply have the same overall parameters as the input sequence.

**Example**

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.new_var(test)
    MEHPTHDQYDQNHKQEPTGSNPNGTPHETNPQP


### The constant_residue_var()

The ``constant_residue_var()`` generates a new sequence variant where specific residues are held constant. The variant will have different residues than the input sequence (with the exception of those specified to be held constant), **but the overall parameters will remain the same between the variant and the original sequence**. You can specify more than one residue to be held constant.


**Example with one residue constant**

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_residue_var(test, constant=['T'])
    GEPNSQEQHDTTNRGDSQHPIQNNQNPDMPSHN

**Example with two residues constant**

    test = 'QEQNGVDQQETTPRQDYPGNQQPNQQAEGQQMQ'
    create.constant_residue_var(test, constant=['T', 'Q'])
    QEQSANDQQETTPKQEAPSPQQASQQHEGQQPQ


### The shuffle_var()

The ``shuffle_var()`` generates a variant that will shuffle specific regions of an IDR. Multiple regions can be specified simultaneously. Regions being shuffled cannot overlap with each otehr
**Note** - The shuffle_var does **NOT** use index values like you would normally in Python. For the shuffle_var, 1 = the first amino acid in the sequence **NOT 0**. 

**Example with one shuffled region**

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.shuffle_var(test, shuffle=[3,9])
    QQNQENNEEDDDQQQEEENNNDDD

**Example with two residues constant**

    test = 'QQQEEENNNDDDQQQEEENNNDDD'
    create.shuffle_var(test, shuffle=[[3,9], [15, 23]])
    QQNNEEQNEDDDQQNDENNEDEQD

**Notice that when you specify 2 regions, you use a list of lists (a nested list).**

**The following variants can change the overall parameters of the generated sequence, unlike previous variants.**


### The constant_class_hydro_var()

The ``constant_class_hydro_var()`` makes a variant that adjusts the hydropathy only via switching between amino acids with different hydropathy values within the same class. Therefore, the order of amino acids *by class* will be the same in the returned variant, but the hydropathy will be adjusted. There is only a limited extent to which the hydropathy van be altered due to the fact that amino acid classes cannot be changed. If you try to change to a value outside of the possible hydropathy values, the error message will include the minimum and maximum possible theoretical values.


**Example decreasing hydropathy** - 
The starting hydropathy of the sequence below is  2.0272. Let's raise it to around 2.7.

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.constant_class_hydro_var(test, hydropathy=2.7)
    GTGGTKIETKTEKKGETTHKTTHTDGLKHTDKKKTHDKSAASRE

**Example where hydropathy is raised higher than possible**

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.constant_class_hydro_var(test, hydropathy=3.7)
    goose.goose_exceptions.GooseInputError:

    Unable to get to objective hydropathy without changing classes of residues.
	For this sequence the lowest possible hydrpathy is 1.611364.
	For this sequence the highest possible hydropathy is 2.834091.


### The kappa_var()

The ``kappa_var()`` variant allows you to change the charge asymmetry of a sequence **provided it has postively and negatively charged residues.** Higher kappa values increase charge asymmetry, lower kappa values reduce charge asymmetry. Values can be between 0 and 1. 

**Example** - 

First we can take something with very symmetrically positions oppositely charged amino acids and increase the kappa value. For reference, the starting kappa value for this 'test' sequence was 0.0012.

    test = 'QNEKRDQNEKRDQNEKRDQNEKRDQNEKRDQN'
    create.kappa_var(test, kappa=0.9)
    KRRRRKKKRKQNQNQNQNQNEEDQDDEDDEEN

Now we can take this newly generated, asymmetrically positioned opposite charged variant and make the charges more moderately symmetrical (something between what we started with and what we made in the previous example).

    previous_variant = 'KRRRRKKKRKQNQNQNQNQNEEDQDDEDDEEN'
    create.kappa_var(previous_variant, kappa=0.35)
    QEEEEDDDRDREKDKKRKNQNQNQNQNQRRKN

**note** GOOSE will allow deviation from your input kappa value by up to 0.02. This is to keep GOOSE from being extremely slow. If you need something closer to your desired value, you can try generating a few variants. You'll likely quickly get what you want.


### The minimal_var()

The ``minimal_var()`` variant allows you to input a sequence of interest and then choose parameters including hydropathy, fcr, ncpr, and scd (sequence charge decoration, which is another measurement of sequence charge asymmetry(SCD)) to alter in your returned sequence. SCD will eventually be changed over to kappa. The objective of this variant is to generate a sequence with the desired input parameters while changing as few amino acids in the sequence as possible. The more you change a given parameter, the more the returned variant sequence will differ from the original.


**Changing one parameter** - 

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.minimal_var(test, hydropathy=3)
    AGAGGRAEGRGERKGEGGGKGGAGDGARGGDRRRGGDKGAAGRE

**Changing multiple parameters** - 

    test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
    create.minimal_var(test, hydropathy=3, fcr=0.2)
    GSGGSRAENRTEQQGEQTSSSSSSGGARSTTRRSSSSKSAASRG


## About sequence analysis in GOOSE

GOOSE provides powerful sequence analysis tools including eight machine learning-based predictors that utilize a LSTM BRNN architecture (nine if you included metapredict v2). All of this functionality is in the 'analyze' module. To use this module simply start by importing it

    from goose import analyze

and then you have a full suite of anlysis tools. The analyze module includes tools for calculating and predicting the following characterstics of your protein of interest (which doesn't even need to be an IDR made in GOOSE!). 

1. Sequence length
2. Fraction of charged residues (FCR)
3. Net charge per residue (NCPR)
4. Average hydropathy
5. Sigma is a property that quantifies charge asymmetry and is defined as the NCPR^2 / FCR
6. delta is a summation of sigma across a sequence using a blob length of 5.5. Basically, more stuff to do with charge asymmetry.
7. kappa - A measurement of charge asymmetry where the max value (1) is the greatest possible charge asymmetry for a given sequence and the min value (0) is the most symmetrical positions possible for oppositely charged residues.
8. Predicted phosphosites for S, T, and Y phosphorylation sites
9. Predicted cellular localization signals including those for nuclear localization, nuclear export and mitochondrial targeting. We will add more as data becomes available.
10. Predicted transcriptional activation domains
11. The fractions of all amino acids in the sequence.

## Using the analyze module in GOOSE

There are 5 functions available in the analyze module after you import it from GOOSE.

### Analyzing general properties

To analyze general properties in GOOSE, you can use the 'properties' function.
This function returns a dict containing all basic properties calculated including length, FCR, NCPR, hydropathy, kappa, and fractions of amino acids.

**Example**

	test = 'GNGGNRAENRTERKGEQTHKSNHNDGARHTDRRRSHDKNAASRE'
	analyze.properties(test)
	{'length': 44, 'FCR': 0.409091, 'NCPR': 0.090909, 'hydropathy': 2.027273, 'kappa': 0.015391434577958246, 'fractions': {'A': 0.09091, 'C': 0.0, 'D': 0.06818, 'E': 0.09091, 'F': 0.0, 'G': 0.11364, 'H': 0.09091, 'I': 0.0, 'K': 0.06818, 'L': 0.0, 'M': 0.0, 'N': 0.13636, 'P': 0.0, 'Q': 0.02273, 'R': 0.18182, 'S': 0.06818, 'T': 0.06818, 'V': 0.0, 'W': 0.0, 'Y': 0.0}}


### Predicting phosphosites

GOOSE has 3 separate networks each trained on a different phosphorylation event (S, T, and Y phosphosites). The phosphosites() function in analyze gives you information on all of them at once. 

**Example**

	test = 'GNGGNRAENRTSSKSERKGEQTHKSNHNDGARHTDRRRSHYDKNAASRE'
	analyze.phosphosites(test)
	{'S': [11, 12], 'Y': [40], 'T': [21, 33]}

In the returned dict, the S is for serine phosphosites, Y for tyrosine phosphosites, and so on. If GOOSE does not predict a phosphosite for a location, it will not return anything for that amino acid. **It is extremely important to note that this predictor is not a gaurentee of a phosphorylation event!**. Protein phospohrylation is incredibly complex, this predictor should be used more as a way to check on something that you want to avoid being phosphorylated (although as with any 'predictor', nothing can be gaurenteed 100%).


### Predicting subcellular localization

If you design an IDR and it ends up somewhere you don't want it to, that's a bad day in the lab. To try to mitigate this problem, we are working on machine learning-based predictors of cellular localization sequences to try to determine where a given protein might end up. As stated in the phosphosite section, this is not a perfect solution, but it is nonetheless better than nothing. This function uses three separate LSTM BRNN networks: 1. NLS, 2. NES, 3. mitochondrial targeting sequences. 

**Example**

For a known mitochondrial protein...

	test = 'MAAAAASLRGVVLGPRGAGLPGARARGLLCSARPGQLPLRTPQAVALSSKSGLSRGRKVMLSALGMLAAGGAGLAMALHS'
	analyze.cellular_localization(test)
	{'mitochondria': {'MAAAAASLRGVVLGPRGAGLPGARARGLLCSARPGQLPLRTPQAVALSSKSGLSRGRKVMLSAL': [1, 65]}, 'NES': 'No NES sequences predicted.', 'NLS': 'No NLS targeting sequences predicted.'}

In the returned dict, the key 'mitochondria' brings up the sequence predicted to be the targeting sequence as well as the coordinates of that sequence (where 1 is the first amino acid). For the other two keys, 'NES' for nuclear export sequences and 'NLS' for nuclear localization sequences, because none were detected the value for each of those key value pairs just states that none were predicted.


### Predicting transcriptional activation domains

If you design an IDR that could bind to DNA, you might not to want it to inadvertantly have a transcriptional activation domain (TAD). You can check that using GOOSE. 

**Example**

For a subset of a protein with a known TAD...

	test = 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'
	analyze.transcriptional_activation(test)
	{'TGNLNSSYFNSLNIDSML': [31, 49]}
	
If there is a TAD present, the function returns the TAD subsequence along with the coordinates for the TAD in the input sequence. Importantly, once again, this is just a predictor and does not gaurentee that a TAD is present or absent from your sequence.


### Predict everything

If you just want a summary of... well basically everything we've covered so far from properties to all predicted features, that's pretty easy to do! 

**Example**

	test = 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'
	analyze.everything(test)
	{'length': 69, 'FCR': 0.202899, 'NCPR': -0.115942, 'hydropathy': 3.413043, 'sigma': 0.066252, 'delta': 0.057173, 'SCD': 0.78989, 'kappa': 0.24406233325921825, 'predicted phosphosites': {'S': [67], 'Y': [18, 37, 50], 'T': [30, 64]}, 'predicted cellular localization': {'mitochondria': 'No mitochondrial targeting sequences predicted.', 'NES': 'No NES sequences predicted.', 'NLS': 'No NLS targeting sequences predicted.'}, 'predicted transcriptional activation': {'TGNLNSSYFNSLNIDSML': [31, 49]}, 'fractions': {'A': 0.0, 'C': 0.0, 'D': 0.14493, 'E': 0.01449, 'F': 0.01449, 'G': 0.01449, 'H': 0.0, 'I': 0.02899, 'K': 0.01449, 'L': 0.14493, 'M': 0.01449, 'N': 0.23188, 'P': 0.01449, 'Q': 0.01449, 'R': 0.02899, 'S': 0.21739, 'T': 0.04348, 'V': 0.01449, 'W': 0.0, 'Y': 0.04348}, 'sequence': 'PNNLNEKLRNQLNSDTNSYSNSISNSNSNSTGNLNSSYFNSLNIDSMLDDYVSSDLLLNDDDDDTNLSR'}

The analyze.everything() function will return a dictionary holding all of the information from sequence properties to predicted phosphosites, cellular localization, and transcriptional activation domains all from one simple function!



### Generating sequences with predicted secondary structure

'''
under construction
'''





## How to cite GOOSE

You can not currently cite GOOSE as we have yet to publish it (hopefully soon!). We would appreciate if you would mention GOOSE in your methods section with a link to the Github page so readers of your paper can understand how you generated the sequences you used.


### Copyright

Copyright (c) 2022, Ryan Emenecker - Holehouse Lab


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
