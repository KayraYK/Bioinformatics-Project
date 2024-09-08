#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdbool.h>

#define SIZE 100
#define NUMPROTEINS 64

bool isBasePair (char neu1, char neu2);

bool isItaDnaSequence (char strand1 [SIZE], char strand2 [SIZE]);

void reverse (char aStrand [SIZE]);

void complementIt (char aStrand [SIZE]);

bool isItPalindrome (char aStrand [SIZE]);

bool isStrandDnaPalindrome (char aStrand [SIZE]);

int howMany (char aStrand [SIZE], char neu);

void dnaToMrna (char aSeq [SIZE], char mRNA [SIZE]);

void translateDnaToMrnaProteins (char aSeq [SIZE]);

int main() {
    char strand1[SIZE];
    char strand2[SIZE];

    // Test isBasePair
    printf("Testing isBasePair function:\n");
    printf("Enter two characters for testing base pairs (e.g., A T): ");
    char bp1, bp2;
    scanf("%c %c", &bp1, &bp2);
    printf("isBasePair result: %d\n", isBasePair(bp1, bp2));

    // Test isItaDnaSequence
    printf("\nTesting isItaDnaSequence function:\n");
    printf("Enter two DNA strands for testing DNA sequence (e.g., ATCG TACG): ");
    scanf("%s %s", strand1, strand2);
    printf("isItaDnaSequence result: %d\n", isItaDnaSequence(strand1, strand2));

    // Test reverse
    printf("\nTesting reverse function:\n");
    printf("Enter a DNA strand for reversing (e.g., ATCG): ");
    scanf("%s", strand1);
    printf("Before Reverse: %s\n", strand1);
    reverse(strand1);
    printf("After Reverse: %s\n", strand1);

    // Test complementIt
    printf("\nTesting complementIt function:\n");
    printf("Enter a DNA strand for complementing (e.g., ATCG): ");
    scanf("%s", strand1);
    printf("Before Complement: %s\n", strand1);
    complementIt(strand1);
    printf("After Complement: %s\n", strand1);

    // Test isItPalindrome
    printf("\nTesting isItPalindrome function:\n");
    printf("Enter a DNA strand for testing palindrome (e.g., ATCG): ");
    scanf("%s", strand1);
    printf("isItPalindrome result: %d\n", isItPalindrome(strand1));

    // Test isStrandDnaPalindrome
    printf("\nTesting isStrandDnaPalindrome function:\n");
    printf("Enter a DNA strand for testing DNA palindrome (e.g., ATCG): ");
    scanf("%s", strand1);
    printf("isStrandDnaPalindrome result: %d\n", isStrandDnaPalindrome(strand1));

    // Test howMany
    printf("\nTesting howMany function:\n");
    printf("Enter a DNA strand and a character to count (e.g., ATCG A): ");
    char ch;
    scanf("%s %c", strand1, &ch);
    printf("howMany result: %d\n", howMany(strand1, ch));

    // Test dnaToMrna
    printf("\nTesting dnaToMrna function:\n");
    printf("Enter a DNA sequence for converting to mRNA (e.g., ATCG): ");
    scanf("%s", strand1);
    char mRNA[SIZE];
    dnaToMrna(strand1, mRNA);
    printf("Original DNA: %s, mRNA: %s\n", strand1, mRNA);

    // Test translateDnaToMrnaProteins
    printf("\nTesting translateDnaToMrnaProteins function:\n");
    printf("Enter a DNA sequence to translate to mRNA and proteins (e.g., ATCG): ");
    scanf("%s", strand1);
    translateDnaToMrnaProteins(strand1);

    return 0;
}