**Search Engine**

Created a commandline search engine based on TREC-Deep-Learning-2020 dataset.

Program.cpp contains the code used to create the search engine.
There are 4 main parts:
1. Create raw postings and a document meta file based on the dataset.
2. Use unix sort to sort the raw postings. (Unix sort is used because it's an I/O efficient sorting method.)
3. Create an inverted index and a lexicon that contains information about the inverted index based on the raw postings.
4. Implement efficient methods for searching within the inverted index and generating search results.

To create a search engine, follow these steps:
- Download the trec dataset from https://microsoft.github.io/msmarco/TREC-Deep-Learning-2020 and put it in a folder called Data. 
- Uncomment the portion of the code in the main function of Program.cpp about generating raw postings and comment the rest of the code in the main function.
- Run unix sort on a linux terminal: sort -k 1,1 -k 2,2n -o ./sortedPostings.txt ./rawPostings.txt
- Uncomment the portion of the code in the main function of Program.cpp about consolidating the raw postings and comment the rest of the code in the main function.

Details of the structure of the inverted index and how it was implemented are written in the pdf
