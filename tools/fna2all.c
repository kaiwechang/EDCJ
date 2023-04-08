#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define true 1
#define false 0
typedef char bool;

typedef struct 
{
    unsigned int markerID;
    int markerFamily;
    char *contigName;
    unsigned int contigOrient;
} sGenomeAllFormat;

typedef struct
{
    unsigned int g1ID;
    int g1;
    unsigned int h1ID;
    int h1;
    unsigned int g2ID;
    int g2;
    unsigned int h2ID;
    int h2;
} PSSPA;

typedef struct
{
    unsigned int g1ID;
    int g1;
    char *g1Contig;
    unsigned int h1ID;
    int h1;
    char *h1Contig;
    unsigned int g2ID;
    int g2;
    char *g2Contig;
    unsigned int h2ID;
    int h2;
    char *h2Contig;
    unsigned int connectCase;
} EPSSPA;

typedef struct
{
    unsigned int g;
    unsigned int h;
    unsigned int weight;
} sAdjMatrix;

typedef struct
{
    unsigned int g;
    unsigned int h;
} sMaxPair;

typedef struct
{
    unsigned int numberOfContig;
    unsigned int maxLengthOfContigName;
    char *contigName;
    unsigned int *contigLength;
} sGenomeFNAInformation;

typedef struct _MARKERINFO
{
    int seqID;
    unsigned char strand;
    unsigned int start;
    unsigned int end;
    unsigned int length;
    struct _MARKERINFO *next;
} sMarkerInformationInBlockFile;

typedef struct
{
    int *marker;
    int *start;
    unsigned int indexCounter;
} sContigMarkerInformation;

//#define SIBELIA_BACTERIAL_M (70)
//#define SIBELIA_PLANT_M     (40)
//#define SIBELIA_M           SIBELIA_BACTERIAL_M
#define SIBELIAZ_M          (200)
#define SIBELIAZ_K          (125)
#define SIBELIAZ_B          (500)

int SIBELIA_BACTERIAL_M = 70;
int SIBELIA_PLANT_M = 40;
int SIBELIA_M = 70;

void ReadSibeliaAll (const char *referenceAllPath, 
                     const char *targetAllPath, 
                     sGenomeAllFormat **referenceAll, 
                     sGenomeAllFormat **targetAll, 
                     unsigned int *numberOfMarkerInReference, 
                     unsigned int *numberOfMarkerInTarget, 
                     unsigned int *maxLengthOfReferenceContigName,
                     unsigned int *maxLengthOfTargetContigName)
{
    FILE *inputFilePointer = NULL;
    char buffer[1024] = {};
    unsigned int maxLengthOfContigName = 0;
    int numberOfMarker = 0;
    int i = 0;

    if ((inputFilePointer = fopen (referenceAllPath, "rb")) == NULL)
    {
        printf ("Error: Cannot read %s\n", referenceAllPath);
        exit (0x0006);
    }
    while (!feof (inputFilePointer))
    {
        fscanf (inputFilePointer, "%*d %*d %s %*d\n", buffer);
        if (strlen (buffer) > maxLengthOfContigName)
        {
            maxLengthOfContigName = strlen (buffer);
        }
        numberOfMarker++;
    }
    fseek (inputFilePointer, 0, SEEK_SET);
    *maxLengthOfReferenceContigName = maxLengthOfContigName;    
    *numberOfMarkerInReference = numberOfMarker;
    *referenceAll = (sGenomeAllFormat *)malloc (sizeof (sGenomeAllFormat) * numberOfMarker);
    for (i = 0; i != numberOfMarker; i++)
    {
        (*referenceAll)[i].contigName = (char *)malloc (maxLengthOfContigName + 1);
        memset ((*referenceAll)[i].contigName, 0, maxLengthOfContigName + 1);
    }
    for (i = 0; !feof (inputFilePointer); i++)
    {
        fscanf (inputFilePointer, "%d %d %s %d\n", 
                                  &((*referenceAll)[i].markerID),
                                  &((*referenceAll)[i].markerFamily),
                                  (*referenceAll)[i].contigName,
                                  &((*referenceAll)[i].contigOrient));
    }
    fclose (inputFilePointer);

    numberOfMarker = 0;
    if ((inputFilePointer = fopen (targetAllPath, "rb")) == NULL)
    {
        printf ("Error: Cannot read %s\n", targetAllPath);
        exit (0x0006);
    }
    while (!feof (inputFilePointer))
    {
        fscanf (inputFilePointer, "%*d %*d %s %*d\n", buffer);
        if (strlen (buffer) > maxLengthOfContigName)
        {
            maxLengthOfContigName = strlen (buffer);
        }
        numberOfMarker++;
    }
    fseek (inputFilePointer, 0, SEEK_SET);
    *maxLengthOfTargetContigName = maxLengthOfContigName;
    *numberOfMarkerInTarget = numberOfMarker;
    *targetAll = (sGenomeAllFormat *)malloc (sizeof (sGenomeAllFormat) * numberOfMarker);
    for (i = 0; i != numberOfMarker; i++)
    {
        (*targetAll)[i].contigName = (char *)malloc (maxLengthOfContigName + 1);
        memset ((*targetAll)[i].contigName, 0, maxLengthOfContigName + 1);
    }
    for (i = 0; !feof (inputFilePointer); i++)
    {
        fscanf (inputFilePointer, "%d %d %s %d\n",
                                  &((*targetAll)[i].markerID),
                                  &((*targetAll)[i].markerFamily),
                                  (*targetAll)[i].contigName,
                                  &((*targetAll)[i].contigOrient));
    }
    fclose (inputFilePointer);
}

void ReadGenomeContigName (const char *filePath, sGenomeFNAInformation **genome, unsigned int *maxLengthOfContigName)
{
    FILE *inputFilePointer = NULL;
    char lengthOfContigName = 0;
    char buffer[1024];
    int contigArrayIndex = 0;

    if ((inputFilePointer = fopen (filePath, "rb")) == NULL)
    {
        printf ("Error: cannot read FNA file (%s).\n", filePath);
        exit (0x0003);
    }

    *genome = (sGenomeFNAInformation *)malloc (sizeof (sGenomeFNAInformation));
    /* Scan the FNA file first to obtain the maximum length among the contig names */
    (*genome)->maxLengthOfContigName = 0;
    (*genome)->numberOfContig = 0;
    while (!feof (inputFilePointer))
    {
        fscanf (inputFilePointer, "%s", buffer);
        if (buffer[0] == '>')
        {
            if (buffer[strlen (buffer) - 1] == '\n')
            {
                lengthOfContigName = strlen (buffer) - 2;
            }
            else
            {
                lengthOfContigName = strlen (buffer) - 1;
            }
            if (lengthOfContigName >= (*genome)->maxLengthOfContigName)
            {
                (*genome)->maxLengthOfContigName = lengthOfContigName;
            }

            ((*genome)->numberOfContig)++;
        }
    }
    *maxLengthOfContigName = (*genome)->maxLengthOfContigName;

    /* Scan the FNA file again to record contig names */
    fseek (inputFilePointer, 0, SEEK_SET);
    (*genome)->contigName = (char *)malloc (((*genome)->maxLengthOfContigName) * ((*genome)->numberOfContig));
    memset ((*genome)->contigName, 0, ((*genome)->maxLengthOfContigName) * ((*genome)->numberOfContig));
    (*genome)->contigLength = (unsigned int *)malloc (sizeof (unsigned int) * ((*genome)->numberOfContig));
    memset ((*genome)->contigLength, 0, sizeof (unsigned int) * ((*genome)->numberOfContig));
    while (!feof (inputFilePointer))
    {
        fscanf (inputFilePointer, "%s", buffer);
        if (buffer[0] == '>')
        {
            if (buffer[strlen (buffer) - 1] == '\n')
            {
                lengthOfContigName = strlen (buffer) - 2;
            }
            else
            {
                lengthOfContigName = strlen (buffer) - 1;
            }
            memcpy (&((*genome)->contigName[contigArrayIndex]), &buffer[1], lengthOfContigName);
            contigArrayIndex += (*genome)->maxLengthOfContigName;
        }
    }
    fclose (inputFilePointer);
}

void SwapTwoElement (int *a, int *b)
{
    int c = *b;

    *b = *a;
    *a = c;
}

void ReadBlockCoordPath (const char * outputPath, 
                         sGenomeFNAInformation **referenceFna, 
                         sGenomeFNAInformation **targetFna,
                         sGenomeAllFormat **referenceMarker,
                         sGenomeAllFormat **targetMarker, 
                         unsigned int *numberOfMarkerInReference,
                         unsigned int *numberOfMarkerInTarget)
{
    char blockCoordPath[1024] = {};
    FILE *inputFilePointer = NULL;
    bool scanFlag;
    char inputChar;
    char buffer[1024] = {};
    unsigned int numberOfGenomeContig;
    unsigned int *numberOfMarkerOnContigInReference = NULL;
    unsigned int *numberOfMarkerOnContigInTarget = NULL;
    unsigned int numberOfContigInReference = 0;
    unsigned int numberOfContigInTarget = 0;
    unsigned int markerNumber = 0;
    unsigned int *newMarkerNumber = NULL;
    bool referenceJoinFlag = false;
    bool targetJoinFlag = false;
    sMarkerInformationInBlockFile *root = NULL;
    sMarkerInformationInBlockFile *current = NULL;
    sMarkerInformationInBlockFile *tempPointer = NULL;
    sContigMarkerInformation *referenceTempContigMarker = NULL;
    sContigMarkerInformation *targetTempContigMarker = NULL;
    unsigned int finalNumberOfMarkerInReference = 0;
    unsigned int finalNumberOfMarkerInTarget = 0;
    bool firstData = false;
    long int markerInformationPosition = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    sprintf (blockCoordPath, "%s//blocks_coords.txt", outputPath);
    if ((inputFilePointer = fopen (blockCoordPath, "rb")) == NULL)
    {
        printf ("Error! Cannot read blocks_coords.txt (generated from Sibelia)\n");
        exit (0x0004);
    }

    numberOfGenomeContig = 0;
    numberOfContigInReference = (*referenceFna)->numberOfContig;
    numberOfContigInTarget = (*targetFna)->numberOfContig;
    scanFlag = true;
    fscanf (inputFilePointer, "%*s %*s %*s\n");
    while (scanFlag)
    {
        fscanf (inputFilePointer, "%c", &inputChar);
        /* read the end of the contig names of reference and target */
        if (inputChar == '-')
        {
            /* read the whole separated symobles in the current line */
            fscanf (inputFilePointer, "%*s\n");
            scanFlag = false;
            continue;
        }
        /* go back one character for reading the correct Seq_id */
        fseek (inputFilePointer, -1, SEEK_CUR);

        if (numberOfGenomeContig < numberOfContigInReference)
        {
            fscanf (inputFilePointer, "%*d %d %s\n", &((*referenceFna)->contigLength[numberOfGenomeContig]), buffer);
            if (buffer[strlen (buffer) - 1] == '\n')
            {
                buffer[strlen (buffer) - 1] = '\0';
            }
        }
        else
        {
            fscanf (inputFilePointer, "%*d %d %s\n", &((*targetFna)->contigLength[numberOfGenomeContig - numberOfContigInReference]), buffer);
            if (buffer[strlen (buffer) - 1] == '\n')
            {
                buffer[strlen (buffer) - 1] = '\0';
            }
        }

        numberOfGenomeContig++;
    }

    numberOfMarkerOnContigInReference = (unsigned int *)malloc (sizeof (unsigned int) * numberOfContigInReference);
    memset (numberOfMarkerOnContigInReference, 0, sizeof (unsigned int) * numberOfContigInReference);
    numberOfMarkerOnContigInTarget = (unsigned int *)malloc (sizeof (unsigned int) * numberOfContigInTarget);
    memset (numberOfMarkerOnContigInTarget, 0, sizeof (unsigned int) * numberOfContigInTarget);

    markerInformationPosition = ftell (inputFilePointer);
    while (!feof (inputFilePointer))
    {
        fscanf (inputFilePointer, "%c", &inputChar);
        /* the marker information */
        if (inputChar == 'B')
        {
            /* obtain the marker number */
            fscanf (inputFilePointer, "%*s %*c%d\n%*s %*s %*s %*s %*s\n", &markerNumber);
            referenceJoinFlag = false;
            targetJoinFlag = false;
            firstData = true;
            scanFlag = true;
            while (scanFlag)
            {
                fscanf (inputFilePointer, "%c", &inputChar);
                if (inputChar == '-')
                {
                    fscanf (inputFilePointer, "%*s");
                    scanFlag = false;
                    continue;
                }
                fseek (inputFilePointer, -1, SEEK_CUR);
                if (firstData)
                {
                    root = (sMarkerInformationInBlockFile *)malloc (sizeof (sMarkerInformationInBlockFile));
                    root->next = NULL;
                    current = root;
                    firstData = false;
                }
                else
                {
                    current->next = (sMarkerInformationInBlockFile *)malloc (sizeof (sMarkerInformationInBlockFile));
                    current = current->next;
                    current->next = NULL;
                }
                fscanf (inputFilePointer, "%d %c %d %d %d\n", 
                                          &(current->seqID), 
                                          &(current->strand), 
                                          &(current->start), 
                                          &(current->end), 
                                          &(current->length));
                if ((current->seqID <= numberOfContigInReference) && (current->length >= SIBELIA_M))
                {
                    referenceJoinFlag = true;
                }
                else
                {
                    targetJoinFlag = true;
                }
            }
            if ((referenceJoinFlag) && (targetJoinFlag))
            {
                while (root != NULL)
                {
                    if (root->length >= SIBELIA_M)
                    {
                        if (root->seqID <= numberOfContigInReference)
                        {
                            numberOfMarkerOnContigInReference[root->seqID - 1]++;
                        }
                        else
                        {
                            numberOfMarkerOnContigInTarget[root->seqID - numberOfContigInReference - 1]++;
                        }
                    }
                    tempPointer = root;
                    root = root->next;
                    free (tempPointer);
                }
            }
            else
            {
                while (root != NULL)
                {
                    tempPointer = root;
                    root = root->next;
                    free (tempPointer);
                }
            }
        }
    }

    newMarkerNumber = (unsigned int *)malloc (sizeof (unsigned int) * (markerNumber + 1));
    memset (newMarkerNumber, 0, sizeof (unsigned int) * (markerNumber + 1));

    referenceTempContigMarker = (sContigMarkerInformation *)malloc (sizeof (sContigMarkerInformation) * numberOfContigInReference);
    for (i = 0; i != numberOfContigInReference; i++)
    {
        referenceTempContigMarker[i].marker = (int *)malloc (sizeof (int) * numberOfMarkerOnContigInReference[i]);
        referenceTempContigMarker[i].start = (int *)malloc (sizeof (int) * numberOfMarkerOnContigInReference[i]);
        referenceTempContigMarker[i].indexCounter = 0;
    }
    
    targetTempContigMarker = (sContigMarkerInformation *)malloc (sizeof (sContigMarkerInformation) * numberOfContigInTarget); 
    for (i = 0; i != numberOfContigInTarget; i++)
    {
        targetTempContigMarker[i].marker = (int *)malloc (sizeof (int) * numberOfMarkerOnContigInTarget[i]);
        targetTempContigMarker[i].start = (int *)malloc (sizeof (int) * numberOfMarkerOnContigInTarget[i]);
        targetTempContigMarker[i].indexCounter = 0;
    }

    fseek (inputFilePointer, markerInformationPosition, SEEK_SET);
    while (!feof (inputFilePointer))
    {
        fscanf (inputFilePointer, "%c", &inputChar);
        /* the marker information */
        if (inputChar == 'B')
        {
            /* obtain the marker number */
            fscanf (inputFilePointer, "%*s %*c%d\n%*s %*s %*s %*s %*s\n", &markerNumber);
            referenceJoinFlag = false;
            targetJoinFlag = false;
            firstData = true;
            scanFlag = true;
            while (scanFlag)
            {
                fscanf (inputFilePointer, "%c", &inputChar);
                if (inputChar == '-')
                {
                    fscanf (inputFilePointer, "%*s");
                    scanFlag = false;
                    continue;
                }
                fseek (inputFilePointer, -1, SEEK_CUR);
                if (firstData)
                {
                    root = (sMarkerInformationInBlockFile *)malloc (sizeof (sMarkerInformationInBlockFile));
                    root->next = NULL;
                    current = root;
                    firstData = false;
                }
                else
                {
                    current->next = (sMarkerInformationInBlockFile *)malloc (sizeof (sMarkerInformationInBlockFile));
                    current = current->next;
                    current->next = NULL;
                }
                fscanf (inputFilePointer, "%d %c %d %d %d\n",
                                          &(current->seqID),
                                          &(current->strand),
                                          &(current->start),
                                          &(current->end),
                                          &(current->length));
                if ((current->seqID <= numberOfContigInReference) && (current->length >= SIBELIA_M))
                {
                    referenceJoinFlag = true;
                }
                else
                {
                    targetJoinFlag = true;
                }
            }
            if ((referenceJoinFlag) && (targetJoinFlag))
            {
                newMarkerNumber[markerNumber] = 1;
                while (root != NULL)
                {
                    if (root->length >= SIBELIA_M)
                    {
                        if (root->seqID <= numberOfContigInReference)
                        {
                            if (root->start < root->end)
                            {
                                referenceTempContigMarker[root->seqID - 1].marker[referenceTempContigMarker[root->seqID - 1].indexCounter] = markerNumber;
                                referenceTempContigMarker[root->seqID - 1].start[referenceTempContigMarker[root->seqID - 1].indexCounter] = root->start;
                            }
                            else
                            {
                                referenceTempContigMarker[root->seqID - 1].marker[referenceTempContigMarker[root->seqID - 1].indexCounter] = -markerNumber;
                                referenceTempContigMarker[root->seqID - 1].start[referenceTempContigMarker[root->seqID - 1].indexCounter] = root->end;
                            }

                            referenceTempContigMarker[root->seqID - 1].indexCounter += 1;
                        }
                        else
                        {
                            i = root->seqID - numberOfContigInReference - 1;
                            if (root->start < root->end)
                            {
                                targetTempContigMarker[i].marker[targetTempContigMarker[i].indexCounter] = markerNumber;
                                targetTempContigMarker[i].start[targetTempContigMarker[i].indexCounter] = root->start;
                            }
                            else
                            {
                                targetTempContigMarker[i].marker[targetTempContigMarker[i].indexCounter] = -markerNumber;
                                targetTempContigMarker[i].start[targetTempContigMarker[i].indexCounter] = root->end;
                            }

                            targetTempContigMarker[i].indexCounter += 1;
                        }
                    }
                    tempPointer = root;
                    root = root->next;
                    free (tempPointer);
                }
            }
            else
            {
                while (root != NULL)
                {
                    tempPointer = root;
                    root = root->next;
                    free (tempPointer);
                }
            }
        }
    }

    fclose (inputFilePointer);

    for (i = 0, j = 1; i <= markerNumber; i++)
    {
        if (newMarkerNumber[i] != 0)
        {
            newMarkerNumber[i] = j++;
        }
    }

    for (i = 0; i != numberOfContigInReference; i++)
    {
        k = referenceTempContigMarker[i].indexCounter;
        finalNumberOfMarkerInReference += k;
        for (j = k - 1; j > 0; j--)
        {
            for (l = 0; l != j; l++)
            {
                if (referenceTempContigMarker[i].start[l] > referenceTempContigMarker[i].start[l + 1])
                {
                    SwapTwoElement (&(referenceTempContigMarker[i].marker[l]), &(referenceTempContigMarker[i].marker[l + 1]));
                    SwapTwoElement (&(referenceTempContigMarker[i].start[l]), &(referenceTempContigMarker[i].start[l + 1]));
                }
            }
        }
    }

    for (i = 0; i != numberOfContigInTarget; i++)
    {
        k = targetTempContigMarker[i].indexCounter;
        finalNumberOfMarkerInTarget += k;
        for (j = k - 1; j > 0; j--)
        {
            for (l = 0; l != j; l++)
            {
                if (targetTempContigMarker[i].start[l] > targetTempContigMarker[i].start[l + 1])
                {
                    SwapTwoElement (&(targetTempContigMarker[i].marker[l]), &(targetTempContigMarker[i].marker[l + 1]));
                    SwapTwoElement (&(targetTempContigMarker[i].start[l]), &(targetTempContigMarker[i].start[l + 1]));
                }
            }
        }
    }

    *referenceMarker = (sGenomeAllFormat *) malloc (sizeof (sGenomeAllFormat) * finalNumberOfMarkerInReference);
    *numberOfMarkerInReference = finalNumberOfMarkerInReference;
    k = 1;
    for (i = 0; i != numberOfContigInReference; i++)
    {
        if (referenceTempContigMarker[i].indexCounter != 0)
        {
            for (j = 0; j != referenceTempContigMarker[i].indexCounter; j++)
            {
                (*referenceMarker)[k - 1].markerID = k;
                (*referenceMarker)[k - 1].contigName = (char *)malloc ((*referenceFna)->maxLengthOfContigName);
                memcpy ((*referenceMarker)[k - 1].contigName, &(*referenceFna)->contigName[i * ((*referenceFna)->maxLengthOfContigName)], (*referenceFna)->maxLengthOfContigName);
                (*referenceMarker)[k - 1].contigOrient = 1;
                if (referenceTempContigMarker[i].marker[j] > 0)
                {
                    (*referenceMarker)[k - 1].markerFamily = newMarkerNumber[referenceTempContigMarker[i].marker[j]];
                }
                else
                {
                    (*referenceMarker)[k - 1].markerFamily = -newMarkerNumber[-(referenceTempContigMarker[i].marker[j])];
                }

                k++;
            }
        }
    }

    *targetMarker = (sGenomeAllFormat *) malloc (sizeof (sGenomeAllFormat) * finalNumberOfMarkerInTarget);
    *numberOfMarkerInTarget = finalNumberOfMarkerInTarget;
    k = 1;
    for (i = 0; i != numberOfContigInTarget; i++)
    {
        if (targetTempContigMarker[i].indexCounter != 0)
        {
            for (j = 0; j != targetTempContigMarker[i].indexCounter; j++)
            {
                (*targetMarker)[k - 1].markerID = k;
                (*targetMarker)[k - 1].contigName = (char *)malloc ((*targetFna)->maxLengthOfContigName);
                memcpy ((*targetMarker)[k - 1].contigName, &(*targetFna)->contigName[i * ((*targetFna)->maxLengthOfContigName)], (*targetFna)->maxLengthOfContigName);
                (*targetMarker)[k - 1].contigOrient = 1;
                if (targetTempContigMarker[i].marker[j] > 0)
                {
                    (*targetMarker)[k - 1].markerFamily = newMarkerNumber[targetTempContigMarker[i].marker[j]];
                }
                else
                {
                    (*targetMarker)[k - 1].markerFamily = -newMarkerNumber[-(targetTempContigMarker[i].marker[j])];
                }

                k++;
            }
        }
    }

    free (numberOfMarkerOnContigInReference);
    free (numberOfMarkerOnContigInTarget);
    free (newMarkerNumber);
    for (i = 0; i != numberOfContigInReference; i++)
    {
        free (referenceTempContigMarker[i].marker);
        free (referenceTempContigMarker[i].start);
    }
    free (referenceTempContigMarker);
    for (i = 0; i != numberOfContigInTarget; i++)
    {
        free (targetTempContigMarker[i].marker);
        free (targetTempContigMarker[i].start);
    }
    free (targetTempContigMarker);
}

void TransformSibeliaToAll (const char * referenceFnaPath, 
                            const char * targetFnaPath, 
                            sGenomeAllFormat **referenceMarkerInformation,
                            sGenomeAllFormat **targetMarkerInformation,
                            unsigned int *numberOfMarkerInReference,
                            unsigned int *numberOfMarkerInTarget,
                            unsigned int *maxLengthOfReferenceContigName, 
                            unsigned int *maxLengthOfTargetContigName, 
                            const char * outputPath)
{
    sGenomeFNAInformation *referenceFna;
    sGenomeFNAInformation *targetFna;
    FILE *outputFilePointer = NULL;
    char path[1024] = {};
    int i = 0;

    //referenceFna = readGenomeContigName (referenceFnaPath);
    ReadGenomeContigName (referenceFnaPath, &referenceFna, maxLengthOfReferenceContigName);
    #ifdef DEBUG
    printf ("Reference has %d contigs and the max length of contigs is %d\n", referenceFna->numberOfContig, referenceFna->maxLengthOfContigName);
    #endif
    ReadGenomeContigName (targetFnaPath, &targetFna, maxLengthOfTargetContigName);
    #ifdef DEBUG
    printf ("Target has %d contigs and the max length of contigs is %d\n", targetFna->numberOfContig, targetFna->maxLengthOfContigName);
    #endif
    ReadBlockCoordPath (outputPath, &referenceFna, &targetFna, referenceMarkerInformation, targetMarkerInformation, numberOfMarkerInReference, numberOfMarkerInTarget);

    /* output the results to reference.all and target.all */
    sprintf (path, "%s//reference.all", outputPath);
    if ((outputFilePointer = fopen (path, "wb")) == NULL)
    {
        printf ("Error! Cannot create %s//reference.all\n", outputPath);
        exit (0x0005);
    }
    for (i = 0; i != *numberOfMarkerInReference; i++)
        fprintf (outputFilePointer, "%d %d %s %d\n", (*referenceMarkerInformation)[i].markerID, (*referenceMarkerInformation)[i].markerFamily, (*referenceMarkerInformation)[i].contigName, (*referenceMarkerInformation)[i].contigOrient);
    fclose (outputFilePointer);

    sprintf (path, "%s//target.all", outputPath);
    if ((outputFilePointer = fopen (path, "wb")) == NULL)
    {
        printf ("Error! Cannot create %s//target.all\n", outputPath);
        exit (0x0005);
    } 
    for (i = 0; i != *numberOfMarkerInTarget; i++)
    {
        fprintf (outputFilePointer, "%d %d %s %d\n", (*targetMarkerInformation)[i].markerID, (*targetMarkerInformation)[i].markerFamily, (*targetMarkerInformation)[i].contigName, (*targetMarkerInformation)[i].contigOrient);
    }
    fclose (outputFilePointer);
}
int main(int argc, char** argv) {
	if (argc < 4) {
		printf("[Error] Usage:\n>>> ./sibeliaToAll <refFnaPath> <tarFnaPath> <outputPath>\n");
		exit(1);
	}
	char* referenceFnaPath = argv[1];
	char* targetFnaPath = argv[2];
	sGenomeAllFormat* referenceMarkerInformation;
	sGenomeAllFormat* targetMarkerInformation;
	unsigned int numberOfMarkerInReference;
	unsigned int numberOfMarkerInTarget;
	unsigned int maxLengthOfReferenceContigName;
	unsigned int maxLengthOfTargetContigName;
	char* outputPath = argv[3];

	char SibeliaCommand[1024] = {};
	sprintf (SibeliaCommand, "Sibelia -k parameter/paraset_bacterial -m 70 %s %s -o %s", referenceFnaPath, targetFnaPath, outputPath);
	system (SibeliaCommand);
	TransformSibeliaToAll(	referenceFnaPath,					targetFnaPath,
							&referenceMarkerInformation,		&targetMarkerInformation,
							&numberOfMarkerInReference,			&numberOfMarkerInTarget,
							&maxLengthOfReferenceContigName,	&maxLengthOfTargetContigName, 
							outputPath);
	return 0;
}
