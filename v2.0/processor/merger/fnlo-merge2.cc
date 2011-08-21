#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>

#include "fnloTable.h"
#include "entry.h"

using namespace std;

int main(int argc, char** argv)
{
	// check parameters
	if (argc < 3)
	{
		cerr << "Usage: fnlo-merge file1.txt [filex.txt]+ result.txt" << endl;
		return 1;
	}

	// loop over arguments
	int nFiles = argc - 1;
	if (access(argv[nFiles], R_OK) == 0)
	{
		cerr << "Error: Output file " << argv[nFiles] << " already exists!" << endl;
		return 1;
	}
	fnloTable resultTable(argv[nFiles]);
	bool validBlockA1 = false, validBlockA2 = false;
	map<Contrib, int> contribCounter;
	map<Contrib, fnloBlockB*> contribMap;

	int nValidTables = 0, Ncontrib = 0, Nmult = 0, Ndata = 0;
	for (int idxFile = 0; idxFile < nFiles - 1; idxFile++)
	{
		string path(argv[idxFile + 1]);
		// File there?
		if (access(path.c_str(), R_OK) == 0)
		{
			fnloTable table(path);
			table.OpenFileRead();

			// Read and check block A1
			if (table.ReadBlockA1() != 0)
			{
				cerr << "Unable to read block A1, skipping " << path << endl;
				continue;
			}
			fnloBlockA1 *blockA1 = table.GetBlockA1();
			if (blockA1->GetItabversion() != 20000)
			{
				cerr << "File does not use table format V20.000, skipping " << path << endl;
				continue;
			}
			if (!validBlockA1)
			{
				*(resultTable.GetBlockA1()) = *blockA1;
				validBlockA1 = true;
			}
			if (!resultTable.GetBlockA1()->IsCompatible(blockA1))
			{
				cerr << "Non compatible A1 blocks found, skipping " << path << endl;
				continue;
			}

			// Read and check block A2
			if (table.ReadBlockA2() != 0)
			{
				cerr << "Unable to read block A2, skipping " << path << endl;
				continue;
			}
			fnloBlockA2 *blockA2 = table.GetBlockA2();
			if (!validBlockA2)
			{
				*(resultTable.GetBlockA2()) = *blockA2;
				validBlockA2 = true;
			}
			if (!resultTable.GetBlockA2()->IsCompatible(blockA2))
			{
				cerr << "Non compatible A2 blocks found, skipping " << path << endl;
				continue;
			}

			int nblocks = blockA1->GetNcontrib() + blockA1->GetNdata();
			for (int i = 0; i < nblocks; i++)
			{
				if (table.ReadBlockB(i) != 0)
				{
					cerr << "Bad format in a B block found, file " << path << endl;
					return 1;
				}
				fnloBlockB *blockB = table.GetBlockB(i);

				Contrib contribution;
				contribution.FromBlock(blockB);
				if (contribMap.find(contribution) == contribMap.end())
				{
					Nmult += contribution.IAddMultFlag;
					Ndata += contribution.IDataFlag;
					Ncontrib += contribution.IAddMultFlag + (contribution.IContrFlag1 > 0 ? 1 : 0);

					fnloBlockB *newBlockB = new fnloBlockB();
					*newBlockB = *blockB;
					newBlockB->BlockA1 = resultTable.GetBlockA1();
					newBlockB->BlockA2 = resultTable.GetBlockA2();
					contribMap[contribution] = newBlockB;
				}
				else
					contribMap[contribution]->Add(blockB);
				++contribCounter[contribution];
			}
			table.DeleteAllBlockB();
			++nValidTables;
		}
		else
			cerr << "Unable to access file, skipping " << path << endl;
	}

	cout << "Found " << nValidTables << " table file(s)." << endl;
	if (nValidTables == 0)
		exit(1);
	cout << "No of types of contributions: " << contribMap.size() << endl;

	// Write result
	cout << "Write merged results to file " << resultTable.GetFilename() << "." << endl;
	resultTable.OpenFileWrite();
	resultTable.GetBlockA1()->SetNcontrib(Ncontrib);
	resultTable.GetBlockA1()->SetNmult(Nmult);
	resultTable.GetBlockA1()->SetNdata(Ndata);
	resultTable.WriteBlockA1();
	resultTable.WriteBlockA2();

	int idx = 0;
	for (map<Contrib, fnloBlockB*>::const_iterator it = contribMap.begin(); it != contribMap.end(); ++it, ++idx)
	{
		cout << " " << contribCounter[it->first] << " file(s) containing ";
		cout << it->first.GetName1() << " " << it->first.GetName2(resultTable.GetBlockA2()->GetILOord());
		resultTable.CreateBlockB(idx, it->second);
		resultTable.WriteBlockB(idx);
		printf("  (%#4.2g events).\n", (double)resultTable.GetBlockB(idx)->GetNevt());
		delete it->second;
	}
	resultTable.CloseFileWrite();

	return 0;
}
