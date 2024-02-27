#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <vector>
#include <unordered_map>
#include <queue>
#include <math.h>

using namespace std;

struct queryEvaluationForDoc
{
    int docID;
    double bm25Score;
    unordered_map<string, int> termToFreq;
    
};

bool operator<(const queryEvaluationForDoc& q1, const queryEvaluationForDoc& q2)
{
    if (q1.bm25Score == q2.bm25Score)
    {
        return q1.docID < q2.docID;
    }
    return q1.bm25Score < q2.bm25Score;
}

//functions
void processString(string& str);
int createRawPosting(string& str, int docID, ofstream& myfile);
int encodeVarByte(vector<int>& numbers, size_t size, char* buffer);
void consolidatePostings(string& sortedPostings, string& invertedIndexFile, string& lexiconFile);
void parseDataSet(string& path, string& rawPostingFile, string& docMetaFile);
void readAndWriteToFiles(fstream& input, ofstream& output, streamoff& inputOffset, int dataSize);
void resetStates(fstream& f, streamoff& offset, int& bufferSize);
vector<int> decodeVarByte(const char* bytes, size_t size);
pair<vector<vector<string>>, double> readDocMeta(string& fileName);
vector<vector<string>> readLexicon(string& file);
int findList(string& term, vector<vector<string>>& lex);
vector<string> getQueryTerms(string query);
int openList(ifstream& invertedIndex, string& term, vector<vector<string>>& lexicon,
    unordered_map<string, char*>& termToInvertedList, unordered_map<string, vector<string>>& termToLexicon,
    unordered_map<string, vector<int>>& termToListMeta);
vector<int> readListMeta(ifstream& invertedIndex, long long invertedListOffset, long metaDataLen);
vector<int> retrieveFreqBlock(int numBlocks, int docsInLastBlock, int blockNumber, vector<int>& listMeta, char* invertedList);
int binarySearch(int key, vector<int>& vec);
void processQuery(vector<vector<string>>& lexicon, string query, bool isDisjunct, pair<vector<vector<string>>, double>& docMeta, 
    ifstream& invertedIndex, ifstream& sourceFile, int numResults);
int getFreq(char* invertedList, int docID, vector<string>& lexicon, vector<int>& listMeta, 
    pair<int, vector<int>>& currentIDBlock, pair<int, vector<int>>& currentFreqBlock);
double computeBM25(int totalNumDoc, int ft, int fdt, int lenDoc, double avgLenDoc);
vector<int> retrieveIDBlock(int numBlocks, int docsInLastBlock, int blockNumber, vector<int>& listMeta, char* invertedList);
int getBlockNumber(int numBlocks, vector<int>& listMeta, int k, int start);
int getNextGEQFromBlock(int k, vector<int>& currBlock);
int nextGEQ(char* invertedList, int k, vector<string>& lexicon, vector<int>& listMeta, pair<int, vector<int>>& currentIDBlock);
string getDisjunctSnippet(int docID, vector<string>& terms, vector<vector<string>>& docMeta, ifstream& source);
string getConjunctSnippet(int docID, vector<string> terms, vector<vector<string>>& docMeta, ifstream& source);
void startSearch(string& docMetaFile, string& lexiconFile, string& invertedIndexFile, string& sourceFile);
void runQueries(vector<vector<string>>& lexicon, pair<vector<vector<string>>, double>& docMeta, ifstream& invertedIndex, ifstream& sourceFile);
vector<queryEvaluationForDoc> getDisjunctQueryResultsRanked(unordered_map<string, char*>& termToInvertedList, unordered_map<string, vector<string>>& termToLexicon,
    unordered_map<string, vector<int>>& termToListMeta, pair<vector<vector<string>>, double>& docMeta, vector<pair<int, string>>& numDocsToTerm, int numResults);
vector<queryEvaluationForDoc> getConjunctQueryResultsRanked(unordered_map<string, char*>& termToInvertedList, unordered_map<string, vector<string>>& termToLexicon,
    unordered_map<string, vector<int>>& termToListMeta, unordered_map<string, pair<int, vector<int>>>& termToCurrentIDBlock,
    unordered_map<string, pair<int, vector<int>>>& termToCurrentFreqBlock, pair<vector<vector<string>>, double>& docMeta,
    vector<pair<int, string>>& numDocsToTerm, int numResults);
vector<string> getQueryTerms(string query);
vector<queryEvaluationForDoc> getTopResults(vector<pair<int, string>>& numDocsToTerm, vector<pair<double, int>>& impactScores, vector<vector<int>>& termFreq, int numResults);
vector<int> processListMeta(int numBlocks, vector<int>& listMeta);
vector<int> decodeIDBlock(char* buffer, int n, int offset, int firstNum);
vector<int> decodeFreqBlock(char* buffer, int n, int offset);

int main()
{
    //Uncomment the part needed
    //generate docMeta and raw postings using sourceFile
    //string sourceFile = "../Data/fulldocs-new.trec";
    //string rawPostingFile = "./rawPostings.txt";
    //string docMetaFile = "./docMeta.txt";
    //parseDataSet(sourceFile, rawPostingFile, docMetaFile);

    //unix sort 
    //operationSort();
    
    //consolidate sorted postings to an inverted index and generate a lexicon file
    //string sortedPostings = "./sortedPostings.txt";
    //string lexiconFile = "./lexicon.txt";
    //string invertedIndexFile = "./invertedIndex.bin";
    //consolidatePostings(sortedPostings, invertedIndexFile, lexiconFile);
    //
    ////start search engine
    string docMetaFile = "./docMeta.txt";
    string lexiconFile = "./lexicon.txt";
    string invertedIndexFile = "./invertedIndex.bin";
    string sourceFile = "./Data/fulldocs-new.trec";
    startSearch(docMetaFile, lexiconFile, invertedIndexFile, sourceFile);
    return 1;
}

/// <summary>
/// Parses data set at the given path and generates raw postings of the form (term, docID)
/// Generates docMeta file which keeps tracks of the docID, url, doc length, and 
/// offset of the text of the doc in source file
/// </summary>
void parseDataSet(string& path, string& rawPostingFile, string& docMetaFile)
{
    int documentID = 0;
    string tempString;
    string url;

    ifstream txtFile(path, ios::binary);
    int temp = 0; //num of docs to parse
    int tokensInDoc = 0;
    int tempTokens = 0;
    bool isURL = false;
    streamoff txtOffset = 0;
    ofstream rawPostings(rawPostingFile, ios::app |ios::binary);
    ofstream docMeta(docMetaFile, ios::app | ios::binary);

    if (!rawPostings.is_open())
    {
        cout << "cannot open myfile";
    }
    if (!docMeta.is_open())
    {
        cout << "cannot open pageFile";
    }
    if (!txtFile.is_open())
    {
        cout << "cannot open txtFile";
    }

    while (getline(txtFile, tempString))
    {
        if (tempString == "<DOC>")
        {
            continue;
        }
        else if (tempString.substr(0, 7) == "<DOCNO>")
        {
            continue;
        }
        else if (tempString == "<TEXT>")
        {
            isURL = true;
            continue;
        }
        else if (isURL)
        {
            url = tempString;
            isURL = false;
            
            txtOffset = txtFile.tellg();
            
            while (getline(txtFile, tempString) && tempString != "</TEXT>")
            {
                processString(tempString);
                tempTokens = createRawPosting(tempString, documentID, rawPostings);
                tokensInDoc = tempTokens + tokensInDoc;
            }
            getline(txtFile, tempString); //tempString = </DOC>
            docMeta << to_string(documentID) + " " + url + " " + to_string(tokensInDoc) + " " + to_string(txtOffset) + "\n";
            documentID++;
            tokensInDoc = 0;
        }
    }

    rawPostings.close();
    docMeta.close();
    txtFile.close();
}

/// <summary>
/// Converts string to lowercase and changes any punctuation to a space character
/// </summary>
void processString(string& str)
{
    for (int i = 0; i < str.size(); i++) {
        if (ispunct(static_cast<unsigned char>(str[i])))
        {
            str[i] = ' ';
            continue;
        }
        str[i] = (char)tolower(str[i]);
    }
}

/// <summary>
/// Creates raw postings of the form (term, docID) from the given string
/// and writes it to the output file 
/// </summary>
/// <returns>number of postings created</returns>
int createRawPosting(string& str, int docID, ofstream& myfile)
{
    int tokens = 0;
    stringstream ss(str);
    string word;
    string s;
    while (ss >> word) {
        //buff.push_back(word + " " + to_string(docID));
        myfile << word + " " + to_string(docID) + "\n";
        tokens++;
    }
    return tokens;
}

// this function does not work on my system which is in windows, but it should work on a linux computer
// I used the same command directly in wsl
/// <summary>
/// Sorts the raw postings into the output file
/// </summary>
/// <returns>1 if sorting is completed with no problem. 0 otherwise</returns>
int operationSort()
{
    string s = "-k 1,1 -k 2,2n -o ./sortedPostings.txt ./rawPostings.txt";
    if (system(s.c_str()))
    {
        return 1;
    }

    else
    {
        cout << "\n" << "Command processor doesn't exists" << endl;
        return 0;
    }
    return 1;
}

/// <summary>
/// Consolidates sortedPostings into an inverted index and generates the lexicon for the 
/// inverted index in the process
/// </summary>
void consolidatePostings(string& sortedPostings, string& invertedIndexFile, string& lexiconFile)
{
    //read in postings line by line 
    string prevLine = "";
    string currLine = "";
    int actualDocID = -1;
    int finalDocID = -1; //for creating posting
    string term = "";
    int frequency = 0;
    int listLength = 0;

    ifstream postings(sortedPostings);
    ofstream invertedIndex(invertedIndexFile, ios::app | ios::binary); //final inverted index
    ofstream lexicon(lexiconFile, ios::app);

    //temp files for storing partial Inverted list
    fstream tempInvertedIndex("./tempTempInvertedIndex.bin", ios::app | ios::in | ios::out | ios::binary);
    fstream lastBlockIDsFile("./lastBlockIDS.bin", ios::in | ios::out | ios::app | ios::binary);
    fstream idBlockSizesFile("./idBlockSizes.bin", ios::in | ios::out | ios::app | ios::binary);
    fstream freqBlockSizesFile("./freqBlockSizes.bin", ios::in | ios::out | ios::app | ios::binary);


    const int bufferSize = 64;
    vector<int> idBuffer;
    vector<int> freqBuffer;
    vector<int> idBlockSizes;
    vector<int> freqBlockSizes;
    vector<int> lastBlockIDs;

    char* idWriteBuffer = new char[10000];
    char* freqWriteBuffer = new char[10000];
    char* metaWriteBuffer1 = new char[10000];
    char* metaWriteBuffer2 = new char[10000];
    char* metaWriteBuffer3 = new char[10000];

    char* readBlock;
    int numBlocks = 0;

    streamoff tempInvertedIndexOffset = 0;
    streamoff lastBlockIDsOffset = 0;
    streamoff idBlockSizesOffset = 0;
    streamoff freqBlockSizesOffset = 0;

    int tempInvertedIndexBufferSize = 0;
    int lastBlockIDsBufferSize = 0;
    int idBlockSizesBufferSize = 0;
    int freqBlockSizesBufferSize = 0;

    streamoff iiOffset = 0;

    if (!postings.is_open())
    {
        cout << "Unable to open postings file.";
    }


    if (!tempInvertedIndex.is_open())
    {
        cout << "Unable to open file.";
    }

    // get the first posting
    if (getline(postings, currLine))
    {
        stringstream ss(currLine); //each line contains term and doc id 
        string word;
        ss >> word;
        term = word;
        ss >> word; //word = docid 
        actualDocID = stoi(word);
        finalDocID = actualDocID;
        listLength = 1;
        prevLine = currLine;
    }

    while (getline(postings, currLine)) {
        if (currLine == prevLine)
        {
            //same term and docID
            frequency++;
        }
        else
        {
            //add the last docID and freq to buffer 
            idBuffer.emplace_back(finalDocID);
            freqBuffer.emplace_back(frequency);

            stringstream ss(currLine); //each line contains term and doc id 
            string word;
            ss >> word;

            bool lastBlock = (word != term);
            int numDocsInLastBlock = idBuffer.size();
            if (numDocsInLastBlock == bufferSize || lastBlock)
            {
                //if buffer is full or it's the last block of this inverted list
                numBlocks = numBlocks + 1;
                int idBlockSize = encodeVarByte(idBuffer, idBuffer.size(), idWriteBuffer);
                int freqBlockSize = encodeVarByte(freqBuffer, freqBuffer.size(), freqWriteBuffer);

                tempInvertedIndex.write(idWriteBuffer, idBlockSize);
                tempInvertedIndex.write(freqWriteBuffer, freqBlockSize);
                idBuffer.clear();
                freqBuffer.clear();

                tempInvertedIndexBufferSize = tempInvertedIndexBufferSize + idBlockSize + freqBlockSize;

                lastBlockIDs.emplace_back(actualDocID);
                idBlockSizes.emplace_back(idBlockSize);
                freqBlockSizes.emplace_back(freqBlockSize);
                
                if (lastBlockIDs.size() == bufferSize || lastBlock)
                {
                    size_t mdSize1 = encodeVarByte(lastBlockIDs, lastBlockIDs.size(), metaWriteBuffer1);
                    lastBlockIDsFile.write(metaWriteBuffer1, mdSize1);
                    size_t mdSize2 = encodeVarByte(idBlockSizes, idBlockSizes.size(), metaWriteBuffer2);
                    idBlockSizesFile.write(metaWriteBuffer2, mdSize2);
                    size_t mdSize3 = encodeVarByte(freqBlockSizes, freqBlockSizes.size(), metaWriteBuffer3);
                    freqBlockSizesFile.write(metaWriteBuffer3, mdSize3);
                    lastBlockIDs.clear();
                    idBlockSizes.clear();
                    freqBlockSizes.clear();

                    lastBlockIDsBufferSize = lastBlockIDsBufferSize + mdSize1;
                    idBlockSizesBufferSize = idBlockSizesBufferSize + mdSize2;
                    freqBlockSizesBufferSize = freqBlockSizesBufferSize + mdSize3;
                }
            }

            if (!lastBlock)
            {
                //inverted list has not ended
                //update new docid and reset frequency
                ss >> word; //word = docid
                int currDocID = stoi(word);
                finalDocID = currDocID - actualDocID; //docID for posting
                actualDocID = currDocID;
                listLength = listLength + 1;

                //frequency starts from 0
                frequency = 0;
                prevLine = currLine;
                continue;
            }

            //update invertedIndex

            //get offset of invertedIndex
            streampos begin, end;
            invertedIndex.seekp(0, ios::beg);
            begin = invertedIndex.tellp();
            invertedIndex.seekp(0, ios::end);
            end = invertedIndex.tellp();
            iiOffset = end - begin;

            //size of metadata 
            int metaSize = lastBlockIDsBufferSize + idBlockSizesBufferSize + freqBlockSizesBufferSize;

            //read and write metadata to invertedindex File 
            readAndWriteToFiles(lastBlockIDsFile, invertedIndex, lastBlockIDsOffset, lastBlockIDsBufferSize);
            readAndWriteToFiles(idBlockSizesFile, invertedIndex, idBlockSizesOffset, idBlockSizesBufferSize);
            readAndWriteToFiles(freqBlockSizesFile, invertedIndex, freqBlockSizesOffset, freqBlockSizesBufferSize);

            //read and write temp inverted index to invertedIndex File
            readAndWriteToFiles(tempInvertedIndex, invertedIndex, tempInvertedIndexOffset, tempInvertedIndexBufferSize);

            //update lexicon
            lexicon << term + " " +
                to_string(numBlocks) + " " +
                to_string(numDocsInLastBlock) + " " +
                to_string(iiOffset) + " " +
                to_string(metaSize) + "\n";

            if (!lexicon.good())
            {
                lexicon.close();
                cout << "error occured when writing to lexicon";
                break;
            }

            //new inverted list
            term = word;
            ss >> word; //word = docid 
            actualDocID = stoi(word);
            finalDocID = actualDocID;
            listLength = 1;

            resetStates(tempInvertedIndex, tempInvertedIndexOffset, tempInvertedIndexBufferSize);
            resetStates(lastBlockIDsFile, lastBlockIDsOffset, lastBlockIDsBufferSize);
            resetStates(idBlockSizesFile, idBlockSizesOffset, idBlockSizesBufferSize);
            resetStates(freqBlockSizesFile, freqBlockSizesOffset, freqBlockSizesBufferSize);

            numBlocks = 0;
            frequency = 0;
            prevLine = currLine;

        }
    }

    //update last posting
    int numDocsInLastBlock = idBuffer.size();

    //add the last docID and freq to buffer 
    idBuffer.emplace_back(finalDocID);
    freqBuffer.emplace_back(frequency);
    numBlocks = numBlocks + 1;

    int idBlockSize = encodeVarByte(idBuffer, idBuffer.size(), idWriteBuffer);
    int freqBlockSize = encodeVarByte(freqBuffer, freqBuffer.size(), freqWriteBuffer);
    tempInvertedIndex.write(idWriteBuffer, idBlockSize);
    tempInvertedIndex.write(freqWriteBuffer, freqBlockSize);
    tempInvertedIndexBufferSize = tempInvertedIndexBufferSize + idBlockSize + freqBlockSize;

    lastBlockIDs.emplace_back(actualDocID);
    idBlockSizes.emplace_back(idBlockSize);
    freqBlockSizes.emplace_back(freqBlockSize);

    size_t mdSize1 = encodeVarByte(lastBlockIDs, lastBlockIDs.size(), metaWriteBuffer1);
    lastBlockIDsFile.write(metaWriteBuffer1, mdSize1);
    size_t mdSize2 = encodeVarByte(idBlockSizes, idBlockSizes.size(), metaWriteBuffer2);
    idBlockSizesFile.write(metaWriteBuffer2, mdSize2);
    size_t mdSize3 = encodeVarByte(freqBlockSizes, freqBlockSizes.size(), metaWriteBuffer3);
    freqBlockSizesFile.write(metaWriteBuffer3, mdSize3);

    lastBlockIDsBufferSize = lastBlockIDsBufferSize + mdSize1;
    idBlockSizesBufferSize = idBlockSizesBufferSize + mdSize2;
    freqBlockSizesBufferSize = freqBlockSizesBufferSize + mdSize3;

    //get offset of invertedIndex
    streampos begin, end;
    invertedIndex.seekp(0, ios::beg);
    begin = invertedIndex.tellp();
    invertedIndex.seekp(0, ios::end);
    end = invertedIndex.tellp();
    iiOffset = end - begin;

    //size of metadata 
    int metaSize = lastBlockIDsBufferSize + idBlockSizesBufferSize + freqBlockSizesBufferSize;

    //read and write metadata to invertedindex File 
    readAndWriteToFiles(lastBlockIDsFile, invertedIndex, lastBlockIDsOffset, lastBlockIDsBufferSize);
    readAndWriteToFiles(idBlockSizesFile, invertedIndex, idBlockSizesOffset, idBlockSizesBufferSize);
    readAndWriteToFiles(freqBlockSizesFile, invertedIndex, freqBlockSizesOffset, freqBlockSizesBufferSize);

    //read and write temp inverted index to invertedIndex File
    readAndWriteToFiles(tempInvertedIndex, invertedIndex, tempInvertedIndexOffset, tempInvertedIndexBufferSize);

    //update lexicon
    lexicon << term + " " +
        to_string(numBlocks) + " " +
        to_string(numDocsInLastBlock) + " " +
        to_string(iiOffset) + " " +
        to_string(metaSize) + "\n";

    postings.close();
    invertedIndex.close();
    tempInvertedIndex.close();
    lastBlockIDsFile.close();
    idBlockSizesFile.close();
    freqBlockSizesFile.close();
    lexicon.close();
}

/// <summary>
/// Encodes the given numbers into varbyte form and save it in the char pointer. 
/// </summary>
/// <returns>The final encoded length</returns>
int encodeVarByte(vector<int>& numbers, size_t size, char* buffer) {
    size_t offset = 0;
    for (size_t i = 0; i < size; i++) {
        int num = numbers[i];
        while (num >= 128) {
            buffer[offset++] = static_cast<unsigned char>(num & 127);
            num >>= 7;
        }
        buffer[offset++] = static_cast<unsigned char>(num + 128);
    }
    return offset;
}


/// <summary>
/// Reads certain bytes off from the input file and writes it to the output file
/// </summary>
void readAndWriteToFiles(fstream& input, ofstream& output, streamoff& inputOffset, int dataSize)
{
    char* readBlock = new char[dataSize];
    input.seekg(inputOffset);
    input.read(readBlock, dataSize);
    output.write(readBlock, dataSize);
}

/// <summary>
/// Sets the fstream to read from the end off the file, the offset to be the current 
/// length of the file, and the bufferSize to be 0
/// </summary>
void resetStates(fstream& f, streamoff& offset, int& bufferSize)
{
    streampos begin, end;
    f.seekg(0, ios::beg);
    begin = f.tellg();
    f.seekg(0, ios::end);
    end = f.tellg();
    f.seekp(0, ios::end);
    offset = end - begin;
    bufferSize = 0;
}

/// <summary>
/// Given a char pointer pointing to a varbyte encoded data,
/// decode certain length of the data 
/// </summary>
/// <returns>vector of integers which are the decoded values </returns>
std::vector<int> decodeVarByte(const char* bytes, size_t size) {
    std::vector<int> numbers;
    int num = 0;
    int shift = 0;
    for (size_t i = 0; i < size; ++i) {
        unsigned char byte = static_cast<unsigned char>(bytes[i]);
        if (byte < 128)
        {
            num = num + (byte << shift);
            shift = shift + 7;
            continue;
        }
        num = num + ((byte - 128) << shift);
        numbers.push_back(num);
        shift = 0; 
        num = 0; // Reset num for the next integer
        
    }
    return numbers;
}

/// <summary>
/// Reads a file about the docMeta of an inverted index and store info of each doc in 
/// a string vector. Also calculates the average length of a doc
/// </summary>
/// <returns>A pair where pair.first is the docMeta stored in a vector 
/// and pair.second is the avg doc length.</returns>
pair<vector<vector<string>>, double> readDocMeta(string& fileName)
{
    vector<vector<string>> vec;

    ifstream file_in(fileName);
    if (!file_in) {
        cout << "error";
    }
    double avgDocLen = 0.0;
    // change this later
    int doc = 1;
    string line;
    while (getline(file_in, line))
    {
        istringstream ss(line);

        vec.push_back({});
        string s;
        for (int i = 0; i < 4; i++)
        {
            ss >> s;
            if (i == 2)
            {
                avgDocLen += (stoi(s) - avgDocLen) / doc;
                ++doc;
            }
            vec.back().push_back(s);
        }
    }
    pair<vector<vector<string>>, double> p;
    p.first = vec;
    p.second = avgDocLen;
    file_in.close();
    return p;
}

/// <summary>
/// Reads a file that contains the lexicon of an inverted index and store the information
/// for each term in a string vector. 
/// </summary>
/// <returns>A vector of string vectors that contains all the information from the lexicon file. </returns>
vector<vector<string>> readLexicon(string& fileName)
{
    vector<vector<string>> vec;

    ifstream file_in(fileName);
    if (!file_in) {
        cout << "error";
    }
    int n = 0;
    string line;
    while (getline(file_in, line))
    {
        istringstream ss(line);
        vec.push_back({});
        string s;
        while (ss >> s)
        {
            vec.back().push_back(s);
        }
        n++;
    }
    file_in.close();
    return vec;
}

/// <summary>
/// Finds where the info regarding a term is located inside a vector that stores the lexicon.
/// </summary>
/// <returns>The index of the term inside the lexicon vector if it can be found. -1 otherwise. </returns>
int findList(string& term, vector<vector<string>>& lex)
{
    int start = 0;
    int end = lex.size() - 1;
    while (start <= end)
    {
        int middle = start + (end - start) / 2;
        if (term == lex[middle][0]) 
        {
            return middle;
        }
        else if (term < lex[middle][0])
        {
            end = middle - 1;
        }
        else {
            start = middle + 1;
        }

    }

    return -1; //term cannot be found in inverted index
     
}

/// <summary>
/// Finds all terms in a query 
/// </summary>
/// <returns>A vector that contains all terms in the query</returns>
vector<string> getQueryTerms(string query)
{
    istringstream ss(query);
    string s;
    vector<string> terms;
    while (ss >> s)
    {
        processString(s);
        terms.push_back(s);
    }
    return terms;
}

/// <summary>
/// Reads and stores the inverted list for a term into memory. 
/// Also decompresses and stores the meta data and lexicon for the inverted list 
/// </summary>
/// <returns>1 if the term exists in the inverted index. -1 otherwise.</returns>
int openList(ifstream& invertedIndex, string& term, vector<vector<string>>& lexicon, 
    unordered_map<string, char*>& termToInvertedList, unordered_map<string, vector<string>>& termToLexicon,
    unordered_map<string, vector<int>>& termToListMeta)
{
    int i = findList(term, lexicon);
    if (i == -1)
    {
        //term is not in lexicon
        return -1;
    }
    vector<string> curr_vec = lexicon.at(i);
    termToLexicon[term] = curr_vec;
    long long offset = stoll(curr_vec.at(3)); //offset of an inverted list is stored at the third index in the lexicon
    long long listLen = 0;
    if (i + 1 >= lexicon.size())
    {
        //term is at the end of inverted index
        //get the listlen by finding the number of bytes from offset to eof
        invertedIndex.seekg(offset, ios::beg);
        streamoff startFrom = invertedIndex.tellg();
        invertedIndex.seekg(0, ios::end);
        streamoff endAt = invertedIndex.tellg();
        listLen = static_cast<long long>(endAt - startFrom);
    } else {
        vector<string> next_vec = lexicon.at(i + 1);
        listLen = stoll(next_vec.at(3)) - stoll(curr_vec.at(3)); 
    }
    long metaDataLen = stol(curr_vec.at(4));
    vector<int> temp;
    temp = readListMeta(invertedIndex, offset, metaDataLen);
    termToListMeta[term] = processListMeta(stoi(curr_vec[1]), temp);
    char* invertedList = new char[listLen];
    invertedIndex.seekg(offset + static_cast<long long>(metaDataLen), ios::beg);
    invertedIndex.read(invertedList, listLen - static_cast<long long>(metaDataLen));
    termToInvertedList[term] = invertedList;

    return 1; 
}

/// <summary>
/// Decompresses the metadata part of an inverted list
/// </summary>
/// <returns>A vector of integers that is the metadata for an inverted lsit</returns>
vector<int> readListMeta(ifstream& invertedIndex, long long invertedListOffset, long metaDataLen)
{
    size_t mDLen = static_cast<size_t>(metaDataLen);
    char* buffer = new char[mDLen];
    invertedIndex.seekg(invertedListOffset, ios::beg);
    invertedIndex.read(buffer, mDLen);
    vector<int> decodedListMeta;
    decodedListMeta = decodeVarByte(buffer, mDLen);
    return decodedListMeta;
}

/// <summary>
/// Creates a list meta vector from the given list meta vector.
/// The new list meta vector will contain the last docID in every docID 
/// block in the inverted list, the offset from the first docID block to every docID block
/// in the inverted list, and the offset from the first docID block to
/// every frequency block in the inverted list. 
/// </summary>
/// <returns>a new list meta vector</returns>
vector<int> processListMeta(int numBlocks, vector<int>& listMeta)
{
    vector<int> newListMeta;
    for (int i = 0; i < numBlocks; i++)
    {
        //last block ids
        newListMeta.push_back(listMeta[i]);
    }
    int offset = 0;
    for (int i = 0; i < numBlocks; i++)
    {
        //docID offsets
        newListMeta.push_back(offset);
        offset = offset + listMeta[numBlocks + i] + listMeta[2 * numBlocks + i];
    }
    
    offset = 0;
    for (int i = 0; i < numBlocks; i++)
    {
        //freq offsets 
        offset = offset + listMeta[numBlocks + i];
        newListMeta.push_back(offset);
        offset = offset + listMeta[2 * numBlocks + i];
    }
    return newListMeta; 
}

/// <summary>
/// Processes a query given the query, its type, and the lexicon, docMeta, invertedIndex and sourceFile
/// </summary>
void processQuery(vector<vector<string>>& lexicon, string query, bool isDisjunct, 
    pair<vector<vector<string>>, double>& docMeta, ifstream& invertedIndex, ifstream& sourceFile, int numResults)
{
    vector<string> queryTerms;
    queryTerms = getQueryTerms(query);
    
    vector<queryEvaluationForDoc> queryResultsRanked;

    unordered_map<string, char*> termToInvertedList;
    unordered_map<string, vector<string>> termToLexicon;
    unordered_map<string, vector<int>> termToListMeta;

    unordered_map<string, pair<int, vector<int>>> termToCurrentIDBlock;
    unordered_map<string, pair<int, vector<int>>> termToCurrentFreqBlock;

    vector<pair<int, string>> numDocsToTerm; //numDocs, term
    
    for (string term: queryTerms)
    {
        int i = openList(invertedIndex, term, lexicon, termToInvertedList, termToLexicon, termToListMeta);
        if (i == -1)
        {
            if (!isDisjunct)
            {
                cout << "I'm sorry. No document matches your query.\n";
                return;
            }
            numDocsToTerm.push_back(make_pair(0, term));
            continue;
        }
        int numDocs = (stoi(termToLexicon[term][1]) - 1) * 64 + stoi(termToLexicon[term][2]);
        numDocsToTerm.push_back(make_pair(numDocs, term));
    }

    ////sort query terms according to their invertedList length
    sort(numDocsToTerm.begin(), numDocsToTerm.end());

    if (isDisjunct)
    {
        queryResultsRanked = getDisjunctQueryResultsRanked(termToInvertedList, termToLexicon, termToListMeta,
            docMeta, numDocsToTerm, numResults);
    }
    else {
        queryResultsRanked = getConjunctQueryResultsRanked(termToInvertedList, termToLexicon, termToListMeta, termToCurrentIDBlock,
            termToCurrentFreqBlock, docMeta, numDocsToTerm, numResults);
    }
    
    //print snippet result to command line 
    if (queryResultsRanked.size() == 0)
    {
        cout << "I'm sorry. No document matches your query.\n";
        return;
    }
    for (auto& r: queryResultsRanked)
    {
        cout << "================================================================================\n";
        cout << docMeta.first[r.docID][1] + "\nlen of doc:" + docMeta.first[r.docID][2] + "\n";
        for (auto& t : queryTerms)
        {
            cout << t + " : " + to_string(r.termToFreq[t]) + " | ";
        }
        cout << "\nbm25 score:" + to_string(r.bm25Score) + "\n";
        if (isDisjunct)
        {
            cout << getDisjunctSnippet(r.docID, queryTerms, docMeta.first, sourceFile) + "\n";
        }
        else
        {
            cout << getConjunctSnippet(r.docID, queryTerms, docMeta.first, sourceFile) + "\n";
        }       
    }

}

/// <summary>
/// Get the conjunctive query results ranked by its BM25 score.
/// The method use to get the results is document-at-a-time query processing 
/// </summary>
/// <returns>a vector containing the final query results in order</returns>
vector<queryEvaluationForDoc> getConjunctQueryResultsRanked(unordered_map<string, char*>& termToInvertedList, 
    unordered_map<string, vector<string>>& termToLexicon, unordered_map<string, vector<int>>& termToListMeta, 
    unordered_map<string, pair<int, vector<int>>>& termToCurrentIDBlock,
    unordered_map<string, pair<int, vector<int>>>& termToCurrentFreqBlock, 
    pair<vector<vector<string>>, double>& docMeta, vector<pair<int, string>>& numDocsToTerm, int numResults)
{ 
    int did = 0;
    int maxDocID = docMeta.first.size() - 1;
    double avgLenOfDoc = docMeta.second;
    vector<pair<double, int>> impactScores;
    vector<vector<int>> termFreq;
   
    for (int i = 0; i < numDocsToTerm.size(); i++)
    {
        vector<int> v(maxDocID + 1, 0);
        termFreq.push_back(v);
    }
    for (auto& p : numDocsToTerm)
    {
        vector<int> temp1;
        vector<int> temp2;
        termToCurrentIDBlock[p.second] = make_pair(-1, temp1);
        termToCurrentFreqBlock[p.second] = make_pair(-1, temp2);
    }
    while (did <= maxDocID)
    {
        string sTerm = numDocsToTerm[0].second;
        did = nextGEQ(termToInvertedList[sTerm], did, termToLexicon[sTerm], termToListMeta[sTerm], termToCurrentIDBlock[sTerm]);
        if (did < 0)
        {
            break;
        }
        bool noMoreIntersection = false;
        int d = -1;
        for (int i = 1; i < numDocsToTerm.size(); i++)
        {
            string temp = numDocsToTerm[i].second;
            d = nextGEQ(termToInvertedList[temp], did, termToLexicon[temp], termToListMeta[temp], termToCurrentIDBlock[temp]);
            if (d < 0)
            {
                noMoreIntersection = true;
                break;
            }
            else if (d > did)
            {
                break;
            }
        }
        if (noMoreIntersection)
        {
            break;
        }
        if (d > did)
        {
            did = d;
        }
        else
        {
            
            int lenOfDoc = stoi(docMeta.first[did][2]);
            vector<int> freqOfTermInDoc;
            vector<int> numDocsContaingTerm;
            double bm25 = 0.0;
            for (int i = 0; i < numDocsToTerm.size(); i++)
            {
                string temp = numDocsToTerm[i].second;
                int freq = getFreq(termToInvertedList[temp], did, termToLexicon[temp], termToListMeta[temp], 
                    termToCurrentIDBlock[temp], termToCurrentFreqBlock[temp]);
                termFreq[i][did] = freq;
                int fdt = freq;
                int ft = numDocsToTerm[i].first;
                bm25 = bm25 + computeBM25(maxDocID + 1, ft, fdt, lenOfDoc, avgLenOfDoc);
            }
            impactScores.push_back(make_pair(bm25, did));
            did++;
        }
    }
    return getTopResults(numDocsToTerm, impactScores, termFreq, numResults);
}

/// <summary>
/// Get the disjunctive query results ranked by its BM25 score.
/// The method use to get the results is term-at-a-time query processing 
/// </summary>
/// <returns>a vector containing the final query results in order</returns>
vector<queryEvaluationForDoc> getDisjunctQueryResultsRanked(unordered_map<string, char*>& termToInvertedList, 
    unordered_map<string, vector<string>>& termToLexicon, unordered_map<string, vector<int>>& termToListMeta, 
    pair<vector<vector<string>>, double>& docMeta, vector<pair<int, string>>& numDocsToTerm, int numResults)
{
    int numDoc = docMeta.first.size();
    double avgLenOfDoc = docMeta.second;
    vector<pair<double, int>> impactScores;
    vector<vector<int>> termFreq;
    
    double k1 = 1.2;
    double b = 0.75;
    
    for (int i = 0; i < numDoc; i++)
    {
        impactScores.push_back(make_pair(0.0, i));
    }
    unordered_map<int, queryEvaluationForDoc> queryResults;
    for (int i = 0; i < numDocsToTerm.size(); i++)
    {
        vector<int> freqForTerm(numDoc, 0);
        termFreq.push_back(freqForTerm);
        string term = numDocsToTerm[i].second;
        int ft = numDocsToTerm[i].first;
        if (ft == 0)
        {
            continue;
        }
        int numBlocks = stoi(termToLexicon[term][1]);
        int docsInLastBlock = stoi(termToLexicon[term][2]);
        
        for (int blockNumber = 1; blockNumber <= numBlocks; blockNumber++)
        {

            vector<int> docIdsInBlock = retrieveIDBlock(numBlocks, docsInLastBlock, blockNumber, termToListMeta[term], termToInvertedList[term]);
            vector<int> freqInBlock = retrieveFreqBlock(numBlocks, docsInLastBlock, blockNumber, termToListMeta[term], termToInvertedList[term]);
            
            for (int j = 0; j < docIdsInBlock.size(); j++)
            {
                int fdt = freqInBlock[j];
                int docID = docIdsInBlock[j];
                termFreq[i][docID] = fdt;
                int lenOfDoc = stoi(docMeta.first[docID][2]);    
                double iScore = computeBM25(numDoc, ft, fdt, lenOfDoc, avgLenOfDoc);
                impactScores[docID].first = impactScores[docID].first + iScore;
            }   
        }
    }

    return getTopResults(numDocsToTerm, impactScores, termFreq, numResults);
}

/// <summary>
/// Get the top numResults documents for a query by sorting all the candidate
/// documents by their impact scores and create a queryEvaluationForDoc 
/// object for each document in the final result. 
/// </summary>
/// <returns>A vector of queryEvaluationForDoc objects which contain the the docID, bm25 score, and query term frequency
/// of each document in the final result</returns>
vector<queryEvaluationForDoc> getTopResults(vector<pair<int, string>>& numDocsToTerm, vector<pair<double, int>>& impactScores, vector<vector<int>>& termFreq, int numResults)
{
    vector<queryEvaluationForDoc> queryResultsRanked;
    if (impactScores.size() == 0)
    {
        return queryResultsRanked;
    }
    priority_queue <pair<double, int>, vector<pair<double, int>>, less<pair<double, int>>> pq(impactScores.begin(), impactScores.end());

    for (int i = 0; i < numResults; i++)
    {
        if (pq.empty())
        {
            break;
        }
        pair<double, int> p = pq.top();
        pq.pop();

        if (p.first == 0)
        {
            break;
        }
        queryEvaluationForDoc qe;
        qe.docID = p.second;
        qe.bm25Score = p.first;
        for (int j = 0; j < numDocsToTerm.size(); j++)
        {
            qe.termToFreq[numDocsToTerm[j].second] = termFreq[j][p.second];
        }
        queryResultsRanked.push_back(qe);
    }
    return queryResultsRanked;
}

/// <summary>
/// Finds the next posting in inverted list with docID >= k
/// </summary>
/// <returns>The docID of the next posting in the inverted list with docID >= k. -1 if none exists.</returns>
int nextGEQ(char* invertedList, int k, vector<string>& lexicon, vector<int>& listMeta, pair<int, vector<int>>& currentIDBlock)
{
    int next = -1;

    if (currentIDBlock.first != -1 && currentIDBlock.second.back() >= k)
    {
        //binary search within the block
        next = getNextGEQFromBlock(k, currentIDBlock.second);
        return next;
    } 

    int numBlocks = stoi(lexicon.at(1));
    int docsInLastBlock = stoi(lexicon.at(2));
    int block = getBlockNumber(numBlocks, listMeta, k, currentIDBlock.first + 1);
    currentIDBlock.first = block;
    if (block == -1)
    {
        return -1;
    }
    currentIDBlock.second = retrieveIDBlock(numBlocks, docsInLastBlock, block, listMeta, invertedList);
    next = getNextGEQFromBlock(k, currentIDBlock.second);
    return next; 
}

/// <summary>
/// Finds the first docID in currBlock >= k
/// </summary>
/// <returns>The first docID in currBlock >= k. -1 if none exists.</returns>
int getNextGEQFromBlock(int k, vector<int>& currBlock)
{
    int start = 0;
    int end = currBlock.size() - 1;
    int sol = -1;
    while (start <= end)
    {
        int middle = start + (end - start) / 2;
        if (currBlock[middle] == k)
        {
            return k;
        }
        if (currBlock[middle] < k)
        {
            start = middle + 1;
        }
        if (currBlock[middle] > k)
        {
            sol = currBlock[middle];
            end = middle - 1;
        }
    }
    return sol;
}

/// <summary>
/// Gets the block number of the id block that contains the first docId >= k.
/// Block number starts from 1
/// </summary>
/// <returns>The block number of the id block that contains the first docID >= k</returns>
int getBlockNumber(int numBlocks, vector<int>& listMeta, int k, int start)
{
    if (start == numBlocks)
    {
        return -1;
    }
    int end = numBlocks - 1;
    int sol = -2;
    while (start <= end)
    {
        int middle = start + (end - start) / 2;
        if (listMeta[middle] < k)
        {
            start = middle + 1;
        }
        else if (listMeta[middle] == k)
        {
            sol = middle;
            break;
        }
        else {
            sol = middle;
            end = middle - 1;
        }
    }
    return sol + 1;
}

/// <summary>
/// Retrieves the decoded version of the id block given the block number and info about the inverted list
/// </summary>
/// <returns>The id block in vector of integers</returns>
vector<int> retrieveIDBlock(int numBlocks, int docsInLastBlock, int blockNumber, 
    vector<int>& listMeta, char* invertedList)
{
    int numDocs = 64;
    if (numBlocks == blockNumber)
    {
        numDocs = docsInLastBlock;
    }
    int offset = listMeta[numBlocks + blockNumber - 1];
    vector<int> decodedBlock;
    int lastIDBlock = 0;
    if (blockNumber != 1)
    {
        lastIDBlock = listMeta[blockNumber - 2];
    }

    decodedBlock = decodeIDBlock(invertedList, numDocs, offset, lastIDBlock);

    return decodedBlock;
}

/// <summary>
/// Decodes a docID block given the inverted list, 
/// the number of docIDs in the block, the position of the block, 
/// and the number to be added to the first docID. 
/// </summary>
/// <returns>The decoded docID block in a vector form.</returns>
vector<int> decodeIDBlock(char* buffer, int n, int offset, int firstNum)
{
    vector<int> numbers;
    int num = 0;
    int shift = 0;
    while (n > 0)
    {
        unsigned char byte = static_cast<unsigned char>(buffer[offset]);
        offset++;
        if (byte < 128)
        {
            num = num + (byte << shift);
            shift = shift + 7;
            continue;
        }
        num = num + ((byte - 128) << shift);
        num = num + firstNum;
        numbers.push_back(num);
        firstNum = num;
        shift = 0;
        num = 0; // Reset num for the next integer
        n--;
    }
    return numbers;
}

/// <summary>
/// Decodes a frequency block given the inverted list, 
/// number of frequencies in the block, and the position of the block,. 
/// </summary>
/// <returns>The decoded frequency block in a vector form.</returns>
vector<int> decodeFreqBlock(char* buffer, int n, int offset)
{
    vector<int> numbers;
    int num = 0;
    int shift = 0;
    while (n > 0)
    {
        unsigned char byte = static_cast<unsigned char>(buffer[offset]);
        offset++;
        if (byte < 128)
        {
            num = num + (byte << shift);
            shift = shift + 7;
            continue;
        }
        num = num + ((byte - 128) << shift) + 1;
        numbers.push_back(num);
        shift = 0;
        num = 0; // Reset num for the next integer
        n--;
    }
    return numbers;
}

/// <summary>
/// Computes BM25 based on the formula
/// </summary>
/// <returns>BM25 score</returns>
double computeBM25(int totalNumDoc, int ft, int fdt, int lenDoc, double avgLenDoc)
{
    double k1 = 1.2;
    double b = 0.75;
    double K = k1 * ((1 - b) + b * lenDoc / avgLenDoc);
    return log(1 + (double(totalNumDoc) - double(ft) + 0.5) / (ft + 0.5)) * (k1 + 1) * fdt / (K + fdt);
}


/// <summary>
/// Finds the frequency of a term in a document given its inverted list and a specific docid in the list
/// </summary>
/// <returns>the frequency of a term in a document</returns>
int getFreq(char* invertedList, int docID, vector<string>& lexicon, vector<int>& listMeta,
    pair<int, vector<int>>& currentIDBlock, pair<int, vector<int>>& currentFreqBlock)
{
    int blockNum = currentIDBlock.first;
    int index = binarySearch(docID, currentIDBlock.second);

    if (currentFreqBlock.first == blockNum)
    {
        return currentFreqBlock.second[index];
    }
    int numBlocks = stoi(lexicon.at(1));
    int docsInLastBlock = stoi(lexicon.at(2));
    vector<int> freqBlock;
    freqBlock = retrieveFreqBlock(numBlocks, docsInLastBlock, blockNum, listMeta, invertedList);
    currentFreqBlock = make_pair(blockNum, freqBlock);
    return freqBlock[index];
}

/// <summary>
/// Finds the index of a key in a vector
/// </summary>
/// <returns>The index if it can be found. -1 otherwise. </returns>
int binarySearch(int key, vector<int>& vec)
{
    int start = 0;
    int end = vec.size() - 1;
    while (start <= end)
    {
        int middle = start + (end - start) / 2;
        if (vec[middle] == key)
        {
            return middle;
        }
        else if (vec[middle] > key)
        {
            end = middle - 1;
        }
        else
        {
            start = middle + 1;
        }
    }
    return -1;
}

/// <summary>
/// Retrieves the decoded version of the freqblock given info about the inverted list 
/// </summary>
/// <returns>the frequency block in vector form </returns>
vector<int> retrieveFreqBlock(int numBlocks, int docsInLastBlock, int blockNumber, 
    vector<int>& listMeta, char* invertedList)
{
    int numDocs = 64;
    if (numBlocks == blockNumber)
    {
        numDocs = docsInLastBlock;
    }
    int offset = listMeta[numBlocks * 2 + blockNumber - 1];

    vector<int> frequencies;
    frequencies = decodeFreqBlock(invertedList, numDocs, offset);
    return frequencies;
}

/// <summary>
/// Generates a snippet for the given document for a disjunct query by referencing the source file.
/// </summary>
/// <returns>The snippet for the given document</returns>
string getDisjunctSnippet(int docID, vector<string>& terms, vector<vector<string>>& docMeta, ifstream& source)
{
    string snippet;
    long long offset = stoll(docMeta[docID][3]); //offset to the start of the doc
    source.seekg(offset, ios::beg);
    string tempStr;
    int lines = 0;
    bool startSnippet = false;
    while (getline(source, tempStr) && lines < 5)
    {
        if (tempStr == "</TEXT>")
        {
            break;
        }
        
        if (startSnippet)
        {
            snippet = snippet +  tempStr + "\n";
            lines++;
        }
        if (!startSnippet) {
            string oriStr = tempStr;
            processString(tempStr);
            stringstream ss(tempStr);
            string s;
            while (ss >> s)
            {
                if (find(terms.begin(), terms.end(), s) != terms.end())
                {
                    startSnippet = true;
                    snippet = snippet + oriStr + "\n";
                    lines++;
                    break;
                }
            }
        }
    }

    return snippet;
}

/// <summary>
/// Generates a snippet for the given document for a conjunct query by referencing the source file.
/// </summary>
/// <returns>The snippet for the given document</returns>
string getConjunctSnippet(int docID, vector<string> terms, vector<vector<string>>& docMeta, ifstream& source)
{
    string snippet;
    long long offset = stoll(docMeta[docID][3]); //offset to the start of the doc
    source.seekg(offset, ios::beg);
    string tempStr;
    int lines = 0;
    int foundTerm = 0;
    if (terms.size() > 4)
    {
        terms.resize(4);
    }
    while (getline(source, tempStr) && lines < 5)
    {
        if (tempStr == "</TEXT>")
        {
            break;
        }
        if (foundTerm != 0)
        {
            snippet = snippet + tempStr + "\n";
            foundTerm--;
            lines++;
            continue;
        }
        if (terms.size() == 0)
        {
            snippet = snippet + tempStr + "\n";
            lines++;
            continue;
        }
        
        processString(tempStr);
        stringstream ss(tempStr);
        string s;
        while (ss >> s)
        {
            auto itr = find(terms.begin(), terms.end(), s);
            if (itr != terms.end())
            {
                terms.erase(itr);
                foundTerm = foundTerm + 1;
                break;
            }
        }

        if (foundTerm != 0)
        {
            if (lines != 0)
            {
                snippet = snippet + "...";
            }
            snippet = snippet + tempStr + "\n";
            lines++;
        }
    }

    return snippet;
}

/// <summary>
/// Accept and process queries according to users' input
/// </summary>
void runQueries(vector<vector<string>>& lexicon, pair<vector<vector<string>>, double>& docMeta, ifstream& invertedIndex, ifstream& sourceFile)
{
    string query;
    string num;
    int disjunct;
    bool isDisjunct = false;
    int status;
    int numResults = 10;
    while (true)
    {
        cout << "Please input your query\n";
        getline(cin, query);
        cout << "Enter 1 for disjunctive results. Enter 0 otherwise.\n";
        getline(cin, num);
        stringstream  disjunctStream(num);
        disjunctStream >> disjunct;
        if (disjunct == 1)
        {
            isDisjunct = true;
        }
        else
        {
            isDisjunct = false;
        }
        if (query == "")
        {
            continue;
        }
        processQuery(lexicon, query, isDisjunct, docMeta, invertedIndex, sourceFile, numResults);   
        
        cout << "Please enter 0 if you would like to continue. Otherwise enter 1.\n";
        getline(cin, num);
        stringstream  continueStream(num);
        continueStream >> status;
        if (status == 1)
        {
            break;
        }
    }
}

void startSearch(string& docMetaFile, string& lexiconFile, string& invertedIndexFile, string& sourceFile)
{
    pair<vector<vector<string>>, double> docMeta;
    docMeta = readDocMeta(docMetaFile);
    vector<vector<string>> lexicon;
    lexicon = readLexicon(lexiconFile);
    ifstream invertedIndex(invertedIndexFile, ios::binary);
    ifstream source(sourceFile, ios::binary);
    runQueries(lexicon, docMeta, invertedIndex, source);
    
    invertedIndex.close();
    source.close();
}
