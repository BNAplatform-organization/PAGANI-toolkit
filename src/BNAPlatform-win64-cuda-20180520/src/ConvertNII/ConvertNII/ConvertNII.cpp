# include <iostream>
# include <fstream>
# include <cmath>
# include <memory.h>
# include <cstring>
# include "Timer.h"
# include "dirent.h" 
//# include "mask_dt.h" 
# include <time.h>
using namespace std;

const int HdrLen = 352;
double ProbCut = 0.5; 


int main(int argc, char * argv[])
{

	if (argc < 3 || argc > 4) 
	{
		cerr<<"Input format: .\\ConvertNII.exe  Directory_containing_files_to_be_converted  niiFile_for_mask  [ |threshold_for_mask]\n"
			<<"For example: .\\ConvertNII.exe  X_deg.nm  mask.nii  [ |0.4]"<<endl;
		system("pause");
		exit(1);	
	}


	DIR *dp;
	struct dirent *dirp;
	if (NULL == (dp = opendir(argv[1])))
	{
		printf("can't open %s", argv[1]);
		exit (1);
	}
	int FileNumber = 0;
	string filenametmp;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);

		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && (filenametmp.substr(filenametmp.find_last_of('.'),3).compare(".nm") == 0) && (filenametmp.size() - filenametmp.find_last_of('.') - 1 == 2))
		{
			FileNumber++;
		}
	}
	cout<<FileNumber<<" files to be processed."<<endl;
	closedir(dp);
	string *filename = new string[FileNumber];
	dp = opendir(argv[1]);
	long long i = 0;
	while((dirp = readdir(dp)) != NULL)
	{
		filenametmp = string(dirp->d_name);
		if (filenametmp.find_last_of('.') == -1)
			continue;
		if(filenametmp.length()>4 && (filenametmp.substr(filenametmp.find_last_of('.'),3).compare(".nm") == 0 ) && filenametmp.size() - filenametmp.find_last_of('.') - 1 == 2)
		{
			filename[i++] = filenametmp;
		}
	}
	// Parse file name

	string mask_file = string(argv[2]);   //Mask file path
	//char * mask_file = argv[2];
	int x_size, y_size, z_size, total_size;
	ifstream fin_mask(mask_file.c_str(), ios_base::binary);
		//fin.open(mask_file, ios::binary);
	if (!fin_mask.good())
	{	cout<<"Can't open\t"<<mask_file.c_str()<<endl;	return 0; }
	
	short hdr[HdrLen / 2];
	fin_mask.read((char*)hdr, HdrLen);
	x_size = hdr[21];
	y_size = hdr[22];
	z_size = hdr[23];
	total_size = x_size * y_size * z_size;
	float *convert_data_f = new float [total_size];	

	
	float *mask = new float [total_size];
	if (hdr[35] == 2) 
	{
		unsigned char * mask_uc = new unsigned char [total_size];
		fin_mask.read((char*)mask_uc, sizeof(unsigned char) * total_size);
		for (int vm = 0; vm<total_size; vm++)
			mask[vm] = (float) mask_uc[vm];
		delete [] mask_uc;
	}
	else if(hdr[35] == 16)
	{
		fin_mask.read((char *)mask, sizeof(float) * total_size);
	}
	else
	{
		cout<<"mask data-type error, Only the type of unsigned char and float can be handled.\n";
		//system("pause");
		return -1;
	}
	fin_mask.close();

	if (argc == 4)
			ProbCut = atof(argv[3]);
	int k = 0 , Numvertice = 0;
	for (k = 0; k < total_size; k++)
		if ((float) mask[k] >= ProbCut)
			Numvertice++;	
	//cout<<Numvertice<<endl;  //debug
	ofstream flog("BNA_time_log", ios::app);

	for (long long i = 0; i < FileNumber; i++)
	{
		
		string b = string(filename[i]);
		string a = string(argv[1]).append("\\").append(filename[i]);  
		cout<<"\nconvert analysis for "<<a.c_str()<<" ..."<<endl;
		ifstream fin(a.c_str(), ios_base::binary);
		if (!fin.good())
		{	cout<<"Can't open\t"<<a.c_str()<<endl;	return 0;}

		const char * input = b.c_str();
		
		string OutNii = a.substr(0, a.find_last_of('.')).append(".nii");
				
		//ifstream fin(input, ios::binary);
		//if (!fin.good())
		//{	cout<<"Can't open\t"<<input<<endl;	system("pause"); return 0;}
		int N = 0;
		fin.read((char*)&N, sizeof(int));
		//cout<<N<<endl; //debug
		if (N!=Numvertice)
		{	cout<<"Mask file doesn't match to the file "<<input<<endl;
			flog<<"Mask file doesn't match to the file "<<input<<endl;
			cout<<"Continue to convert the next file"<<endl;			
			continue;
		}
		float * fraw_data = new float [N];

		fin.read((char*)fraw_data, sizeof(float) * N);
		fin.close();
				

		//int * iraw_data = (int*) fraw_data;
		//int * convert_data = new int[total_size];
			
		
		int l = 0;
	    if (input[strlen(input) - 2] == 'n')
		{
			hdr[35] = 16;
			hdr[36] = 32;
			for (k = 0; k < total_size; k++)
				if ((float) mask[k] >= ProbCut)
					convert_data_f[k] = fraw_data[l++];
				else
					convert_data_f[k] = 0;
		}
		else 
		{
			cout<<"Unknown file type:\t"<<input<<endl;
			return 0;
		}
		ofstream fout(OutNii.c_str(), ios::binary | ios::out);
		fout.write((char*)hdr, HdrLen);
		fout.write((char*)convert_data_f, sizeof(float) * total_size);
		fout.close();		
		delete []fraw_data;
		
		//delete []mask;
		cout<<"Convert to .nii completed!"<<endl;
	}
	delete []mask;
	delete []convert_data_f;
	flog.close();
	cout<<"==========================================================="<<endl;
	system("pause");
	return 0;
}
