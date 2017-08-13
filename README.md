# ms2geno -- convert ms output to genotypes in files

To get and compile this you can do like so:

```sh
# first clone it from GitHub
git clone https://github.com/eriqande/ms2geno

# then get the submodules
cd ms2geno
git submodule init
git submodule update

# then compile it.  The result is the executable with the 
# operating system appended. For example on my Mac I get
# the binary ms2geno-Darwin
./Compile_ms2geno.sh 


# You can do this to get the full built-in help
ms2geno-Darwin --help-full


# and, if you have Dick Hudson's program ms on your system, you can do a small run like this
ms 800 20  -t 8 -I 2 400 400 75 | ms2geno-Darwin -l 10 -b 150 150 -m 50 50 -u .15 2

```

Of course, if you want to just use a precompiled binary in the repo, 
you can do that too.  



The program description and options that you get from 
doing `ms2geno --help-full` are listed below:

```
ms2geno  --  convert ms output to genotypes in files

Author(s):
	Eric C. Anderson (eric.anderson@noaa.gov)

About the Program:
    This program is designed to convert output from the makesamples (also called ms)
    program, written by Richard Hudson, into data sets of microsatellite or SNP
    genotypes. ms yields sequences of 0s and 1s along chromosomes from the infinite
    sites model---a 1 means the sampled chromosome inherited a certain mutated
    nucleotide and a 0 means that received the non-mutated form of the nucleotide.
    It is straightforward to use such strings of 0s and 1s to simulate
    microsatellites via a stepwise mutation process, or as the background for
    ascertaining SNPs, especially when each locus is assumed unlinked and simulated
    on a separate, independent, coalescent tree. 
    
    This version of ms2geno assumes that each locus is independently segregating. A
    future version could be made which used the recombination simulation features
    of ms to allow for linked markers. However, it does not do that yet! 
    
    The code for reading ms output into simple data structures was taken from one
    of the little ms output analysis programs that Dick Hudson distributes with
    makesamples. Thanks for providing source code Dick! 
    
    Here is how the program works. First you have to simulate from ms. You have to
    make sure that you have twice as many chromosomes (set with the nsam option in
    ms) as you want to have diploid genotypes in your ms2geno output data sets. If
    you run ms with samples from multiple subpopulations using the -I option,
    ms2geno will be able to detect that and partition the individuals appropriately
    into the subpopulations. The howmany option in ms sets how many independent
    coalescent trees are simulated by ms. Each one corresponds to a separate locus.
    If you set howmany=100 in ms, say, but choose to simulate only 5 loci using the
    -l/--num-loci option with ms2geno, then ms2geno will output 100/5=20 separate
    data sets, each with five loci. 
    
    The output goes into files: BaseFile_XX.txt and MixFile_XX.txt where XX is the
    data set number. The format is a sort of quasi-GenePop format that works with
    Steven Kalinowski's program GMA, and with my program gsi_sim. It might even
    work with GenePop. BaseFile and MixFile are just two separate files that 
    
    You can control the genetic variability by using the -t (theta) option to ms.
    With microsatellites, increasing the multi-step mutation rate (see option
    -u/---msat) in ms2geno will also increase the variability. If you are
    simulating SNPs with ms2geno it is a good idea to simulate more than enough
    trees in ms for all the loci because SNPs might not get ascertained from some
    of the simulated trees. You can use the -s option to ms to specify how many
    segregating sites should be available for SNP discovery. The details of SNP
    discovery will be covered elsewhere. 
    
    If the file ms2geno_seeds is present in the current working directory then its
    contents are used to seed the random number generator. If not, then seeds are
    generated from the current time. The current state of the random number
    generator is written to ms2geno_seeds at the end of execution so it may be used
    over and over again in a sane fashion with respect to the random number
    generator. 
    
    This has been an incomplete description, but I just wanted to get something in
    there for guiLiner.

In the following:
	"J" refers to an integer argument to an option

	"R" refers to a real number argument to an option

	"S" refers to a string argument to an option

	"F" refers to a file path argument to an option. For example,
		"datfile.txt" if the file is in the current working directory, or
		something like "~/eriq/Documents/data/datfile.txt" if you want to
		provide a complete file path.  (Beware of spaces in file paths!)

	"D" refers to a directory path argument to an option. For example,
		"data_direcory/" if the directory is in the current working directory, or
		something like "~/eriq/Documents/data_directory/" if you want to
		provide a complete directory path.  Note that the trailing slash should be
		optional, but currently is not.  (ERIC ADD MACROS FOR GETTING FILES AND DIRECTORIES

	"G" refers to a string that gives a (possibly) discontinous range of
		nonnegative integers.  For example:  "1-5,7,9,10-15" specifies
		the integers 1 through 5, 7, 9, and 10 through 15.  There can be no
		whitespace in the string specifying the range, and the numbers must all
		be increasing.  Also, the string cannot start nor end with a comma or a dash.
		Finally, you should not use "-" to denote two ranges without placing any
		commas in between.

	"C" refers to a "constrained" string argument to an option,
		i.e., the argument is a string that may only be drawn from a small
		set of alternatives, as specified in the help-full description.

   ****  Program-description and documentation  parameters  ****

     --help              
        this returns a short list of all program options and associated
        arguments

     --help-full         
        this returns a full list of all program options and associated
        arguments

     --help-nroff        
        this returns a full list of all program options and associated
        arguments using the formatting styles in nroff that give you the look
        of a man page. View the formatted ouput by doing: 'prog --help-nroff
        | nroff -man | more' where prog is the name of the program.

     --help-xml          
        This returns a list of all options in a simple XML format which is
        suitable for input to the guiLiner front end.

     --version           
        prints program version information

     --version-history   
        prints history of different versions

     --command-file      F
        By using this option you can store command line arguments in the file
        named in F. You may have any kind of white space (including line
        endings) in the file. The line endings are treated as just another
        space. Comments may be included in the file by enclosing them within
        a pair of ampersands (the & character). Note that you must have a &
        at the beginning and at the end of the comment. You cannot use just a
        single & to comment to the end of a line. Your comments may spread
        over several lines---they will still be stripped from the resulting
        command line so long as the are enclosed in ampersands. This feature
        is helpful if you have long and complex command lines that you wish
        to store if it makes it easier to read the command line by breaking
        it across multiple lines or if you have certain clusters of options
        that you like to store together as a module. This option may be used
        up to 10000 times. Optional.


   ****  Command Line Switches To Input Options to ms2geno  ****

     --ploidy            J
        Each individual is imagined to carry J copies of each gene. This
        option is not required. The default is 2.

-l , --num-loci          J
        This will be the number of loci included in each data set

-b , --baselines         J1 ... JK
        Sets the number of individuals to be sampled from the K simulated
        populations into the baseline samples. K must be equal to the number
        of populations given using the -I argument to ms. If this option and
        -m/--mixtures are used in the same command line they must have the
        same number of arguments K.

-m , --mixtures          J1 ... JK
        Sets the number of individuals to be sampled from the K simulated
        populations into the mixture samples. K must be equal to the number
        of populations given using the -I argument to ms. If this option and
        -b/--baselines are used in the same command line they must have the
        same number of arguments K.

-o , --write-012         
        If this option is given, then the program will write data in 012
        format. This is only compatible with SNPs and ploidy = 2. 012 format
        produces three files for each Basefile and each Mixfile: a .012 file,
        a .012.indv file and a .012.pos file. This option is not required.

-u , --msat              R1 R2
        R1 is the probability of multistep mutations and R2 is the mean of the
        Poisson distribution that determines the step size of multistep
        mutation. The length of multistep mutations is 2 + a draw from a
        Poisson distribution with mean R2. Applying this option invokes the
        microsatellite/stepwise mutation model. This is the default. The
        default value for R1 is 0 in which case R2 is irrelevant. This
        options cannot be used in conjunction with the -s/--snp option.

-s , --snp-asc           J1 J2 J3
        This says that a SNP will be ascertained if it there are J3 or more
        copies of the minor allele among J2 individuals from population J1.
        This option may be given any number of times; however it should not
        be given multiple times for any single population. You may wish to
        count the number of 0s and 1s in a number of different populations
        and then ascertain on the basis of the number of copies of the minor
        allele across all those populations. This can be achieved by calling
        this option for each of the populations that you want to have
        involved; but for each one use J3=-1. Then give the overall minimum
        number of minor alleles needed using the --all-pops-asc option. We
        also have a quick hack, hijacking this option, that allows us to do
        ascertainment on the basis of finding at least one of every possible
        diploid genotype. In other words, you will only ascertain the SNPs if
        the 00 and 11 homozygotes appear in the ascertainment sample, along
        with the 01 or 10 heterozygotes. To use this option, you let the J3
        argument be -2. Finally, if you want to do the pseudo-AFLP thing then
        you specify the population with the J1 argument and you set the J3
        argument to -3. In this case the ascertainment is not done population
        by population. With the pseudo-AFLP thang, the ascertainment always
        depends on the numbers of phenotypes amongst all the populations
        having J3=-3.

     --all-pops-asc      J
        This option can be used for SNP ascertainment. All the gene copies in
        ascertainment sets specified by the -s/--snp-asc option go into a
        pool and a SNP will be ascertained if the least frequent site allele
        occurs in J or more copies in that pool. Note that this criterion is
        applied for each site after all of the population-specific criteria
        given in the -s/--snp-asc options. To ensure that ascertainment is
        based only on the total number of genes amongst all the ascertainment
        populations you should make sure that the J3 argument to all the
        -s/--snp-asc option calls is -1.

     --all-pops-geno-asc 
        This option can be used when any of the populations have ascertainment
        samples that are specified with using the -s or --snp-asc option with
        the J3 argument equal to -2. Such an option causes a SNP to be
        ascertained in such a population if all three genotypes are observed
        in the population. If you issue this --all-pops-geno-asc option, then
        even if all three genotypes are not seen in any of the ascertainment
        samples, a SNP is still ascertained so long as all three genotypes
        appear amongst the inidividuals in the ascertainment samples that had
        the J3 argument equal to -2.

     --all-pops-pseudo-aflp J1 J2
        This is a bit of a hack. I am busy doing some simulations and I need
        this functionality: This will produce SNPs and it will also produce a
        list of which SNP allele is the recessive allele in a dominant marker
        system like ALPs. With that information, you can then hybridize
        individuals using their SNP genotypes and process those results into
        AFLP phenotypes later. Using this option, a SNP will be ascertained
        if at least J1 recessive (non-band-producing) phenotypes and J2
        band-producing phenotypes are seen from amongst all the individuals
        in the populations that get the -s or --snp-asc option having the J3
        argument equal to -3.

     --no-perm           
        Invoking this option inhibits the default behavior of permuting the
        ascertainment population labels for each SNP. Permuting is the
        default behavior because it seems that you could have some biases if
        you scan through the ascertainment populations in the same order
        every time.

     --ms-pop-override   J1 ... JK
        This odd little option is used to specify to ms2geno the number of
        chromosomes in each population simulated by ms WHEN YOU WANT THAT
        NUMBER TO BE DIFFERENT FROM WHAT WAS SPECIFIED ON THE MS COMMAND
        LINE. Please note that this should be in terms of the number of
        chromosomes and not in terms of the number of diploid individuals.
        This is not an option you normally need to use. The only time it
        should be invoked is when you have simulated a large number of
        chromosomes in a SINGLE population with ms but want you ms2geno to
        carve up the individuals from that single population into separate
        populations. Obviously the number of populations should match what is
        given for baseline and mixture files. ms2geno will check to make sure
        that enough chromosomes were simulated in ms to cover the ones
        requested here.

-f , --ms-file           F
        Should you prefer to have redirected the ms output into a file and
        saved it, you can use that saved file as the input to ms2geno by
        either redirecting it to standard input (with < file.txt) or you can
        name the file on the command line using this option. This is included
        because some operating systems may not have good file redirection
        facilities, and because this makes it much easier to use ms2geno with
        the guiLiner front end.

```
