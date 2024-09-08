/***********************************************************************************
Compiling the program
The program should be compiled using the following flags: -std=c99 -Wall

compiling:
gcc -std=c99 -Wall Functions.c main.c

Running: ./a.out
OR
gcc -std=c99 -Wall Functions.c main.c -o Bioinformatics_Functions_Test
Running the Program: ./assn2
***********************************************************************************/ 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define SIZE 100
#define NUMPROTEINS 64

bool isBasePair (char neu1, char neu2) {

    if ( (neu1 == 'T' && neu2 == 'A') || (neu1 == 'A' && neu2 == 'T') || (neu1 == 'C' && neu2 == 'G') || (neu1 == 'G' && neu2 == 'C') ) { //checks if strand one and strand two have pairs throughout both strands
        return true; // true means all pairs are correct
    } else {
        return false; // false means at least one pair is wrong
    }
}

bool isItaDnaSequence(char strand1[SIZE], char strand2[SIZE]) {

    for (int count = 0; count < SIZE; count++) {
        if ((strand1[count] == 'A' && strand2[count] != 'T') || (strand1[count] == 'T' && strand2[count] != 'A') || (strand1[count] == 'C' && strand2[count] != 'G') || (strand1[count] == 'G' && strand2[count] != 'C')) { // checks if both strands of the dna form a correct full dna sequence
            return false;
        }
    }
    return true;
}

void reverse(char aStrand[SIZE]) {

    int len;

    len = strlen(aStrand); // gets length of string

    for (int count = 0; count < len / 2; count++) { // using the length of the string, sets each value of aStrand to its reverse
        char temp = aStrand[count];
        aStrand[count] = aStrand[len - count - 1];
        aStrand[len - count - 1] = temp;
    }
}

void complementIt (char aStrand [SIZE]) { // goes through each element and complements it using many if statements
    for (int count = 0; count < SIZE; count++) {
        if (aStrand[count] == 'C'){
            aStrand[count] = 'G';
        } else if (aStrand[count] == 'G'){
            aStrand[count] = 'C';
        } else if (aStrand[count] == 'T'){
            aStrand[count] = 'A';
        } else if (aStrand[count] == 'A'){
            aStrand[count] = 'T';
        } 
    }
}

bool isItPalindrome(char aStrand[SIZE]) { // Checks if the array is a palindrome
    
    int len; 

    len = strlen(aStrand);

    char aStrand2[SIZE]; 

    for (int count = 0; count < len; count++) {
        aStrand2[count] = aStrand[count];
    }
    
    aStrand2[len] = '\0'; 

    reverse(aStrand2);

    for (int count2 = 0; count2 < len; count2++) {
        if (aStrand[count2] != aStrand2[count2]) {
            return false;
        }
    }

    return true;
}

bool isStrandDnaPalindrome(char aStrand[SIZE]) { // Checks if the array represents a DNA palindrome
 
    int len = strlen(aStrand);

    char aStrand2[SIZE];

    for (int count = 0; count < len; count++) {
        aStrand2[count] = aStrand[count];
    }

    aStrand2[len] = '\0';

    complementIt(aStrand2);
    reverse(aStrand2);

    for (int count2 = 0; count2 < len; count2++) {
        if (aStrand[count2] != aStrand2[count2]) {
            return false;
        }
    }

    return true;
}

int howMany (char aStrand [SIZE], char neu) { // checks how many of the selected neu is in the array
    
    int characterCount = 0;

    for (int count= 0; count <50; count++){
        if (aStrand [count] == neu) {
            characterCount = characterCount + 1;
        }
    } 
    
    return characterCount;
}

void dnaToMrna (char aSeq [SIZE], char mRNA [SIZE]) { // converts the dna sequence into an mRna sequence

    int len = strlen(aSeq);

    char aSeq2 [SIZE];

    for (int count = 0; count < len; count++) {
        aSeq2[count] = aSeq[count];
    }

    aSeq2[len] = '\0';

    complementIt(aSeq2);

    for (int count2 = 0; count2 < len; count2++) {
        if (aSeq2[count2] == 'T') {
            mRNA[count2] = 'U';
        }
        else {
            mRNA[count2] = aSeq2[count2];
        }
    }
    
    mRNA[len] = '\0';
}

char getCode (char protein [3]) {

    // array allProteins has a list of all 3-triplet amino acids required for this assignment
    
    char allProteins [NUMPROTEINS][SIZE] = {"GCA", "GCC", "GCG", "GCU", "AGA","AGG", "CGA", "CGC","CGG","CGU","GAC", "GAU","AAC","AAU","UGC","UGU","GAA","GAG","CAA","CAG", "GGA", "GGC","GGG","GGU","CAC", "CAU","AUA","AUC","AUU", "UUA", "UUG","CUA","CUC", "CUG", "CUU", "AAA", "AAG","AUG", "UUC","UUU", "CCA", "CCC", "CCG", "CCU", "AGC","AGU","UCA","UCC","UCG", "UCU","ACA","ACC","ACG", "ACU", "UGG","UAC","UAU", "GUA","GUC","GUG", "GUU"};
    
    // array allCodes stores the single-letter code of each triplet given in the above array
    
    char allCodes [NUMPROTEINS] = {'A','A','A','A','R','R','R','R','R','R','D','D','N','N','C','C','E','E','Q','Q', 'G','G','G','G','H','H','I','I','I', 'L','L','L','L','L','L','K','K','M','F','F','P','P','P','P','S','S','S','S','S','S','T','T','T','T', 'W','Y','Y','V','V','V','V'};
    
    for (int i = 0; i < NUMPROTEINS; i++) {
        
       if  (strncmp (protein, allProteins[i], 3) == 0){
            return allCodes [i];
        }
    }
    
    return 'Z';   // to indicate an incorrect code - code that doesn't exist in array allCodes
}

void translateDnaToMrnaProteins(char aSeq[SIZE]) { // takes the dna sequence and uses previous functions and new code to translate to mRna proteins

    char mRNA[SIZE];

    dnaToMrna(aSeq, mRNA);

    printf("DNA: %s\n", aSeq);

    printf("mRNA: %s, which translates to:\n", mRNA);

    int len = strlen(mRNA);

    for (int i = 0; i < len; i += 3) {
        if (i + 2 < len) {
            char temp[4] = {mRNA[i], mRNA[i + 1], mRNA[i + 2], '\0'};

            char aminoAcid = getCode(temp);

            printf("%s : %c\n", temp, aminoAcid);
        }
    }
}