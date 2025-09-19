
#include <fstream>
#include <array>
#include <string>
#include <time.h>
#include <direct.h>		
#include <algorithm>

#include <tuple>

using std::array; using std::remove;
using namespace std;

std::string path = "E:/�츮��_2021����/Test����/1)1st_Image/1st_Image";
std::string path_Save = "E:/�츮��_2021����/Code_c++/2)Noise_reduction";

float mean(float* array, int size) {
	float sum = 0.0;
	for (int i = 0; i < size; i++)
		sum += array[i];
	return sum / (float)size;
}

float min_ws(float* array, int size) {
	float min = std::numeric_limits<float>::max();

	for (int i = 0; i < size; i++)
		if (min>array[i]) {
			min = array[i];
		}
	return min;
}

unsigned short max_ws(unsigned short* array, int size) {
	float max = std::numeric_limits<unsigned short>::min();

	for (int i = 0; i < size; i++)
		if (max<array[i]) {
			max = array[i];
		}
	return max;
}

typedef struct Values {
	unsigned short max;
	unsigned short sensor;
	int padded_size;
	int sizeX_pad;
	int sizeY_pad;
}Values;

unsigned int mod_fuction(unsigned int a, unsigned int q)
{
	unsigned int b = a / q;
	return a - b*q;
}

void writeImageData_slice(unsigned short* processed, int sizeX, int sizeY, string f_name, int ith) {

	int	  String_Size = 200;
	char *File_Name = new char[String_Size];
	char *String_Ouput = new char[String_Size];
	unsigned short* save_slice = new unsigned short[sizeY*sizeX];
	const char* path_s = path.c_str();
	const char* path_Save_s = path_Save.c_str();

	// Writing Voxel Slices
	for (int iny = 0; iny < sizeY; iny++) {
		for (int inx = 0; inx < sizeX; inx++) {
			save_slice[iny *sizeX + inx] = (unsigned short)((processed[iny *sizeX + inx]));
		}
	}
	sprintf(File_Name, "/%04d_%s.raw", ith, f_name);
	strcpy(String_Ouput, path_Save_s);
	strcat(String_Ouput, "/Result");
	_mkdir(String_Ouput);
	strcat(String_Ouput, File_Name);

	FILE *Foutput;
	if ((Foutput = fopen(String_Ouput, "wb")) != NULL)
	{
		fwrite(save_slice, sizeof(unsigned short), sizeX*sizeY, Foutput);
	}
	fclose(Foutput);

	printf("Data Save: %s \n", File_Name);

	delete[] save_slice;
	delete[] File_Name;
	delete[] String_Ouput;
}

void writeImageData_slice_FloattoShort(float* processed, int sizeX, int sizeY, string f_name, int ith) {

	int	  String_Size = 200;
	char *File_Name = new char[String_Size];
	char *String_Ouput = new char[String_Size];
	unsigned short* save_slice = new unsigned short[sizeY*sizeX];
	const char* path_s = path.c_str();
	const char* path_Save_s = path_Save.c_str();

	// Writing Voxel Slices
	for (int iny = 0; iny < sizeY; iny++) {
		for (int inx = 0; inx < sizeX; inx++) {
			save_slice[iny *sizeX + inx] = (unsigned short)((processed[iny *sizeX + inx]));
		}
	}

	sprintf(File_Name, "/%04d_%s.raw", ith, f_name);
	strcpy(String_Ouput, path_Save_s);
	strcat(String_Ouput, "/Result");
	_mkdir(String_Ouput);
	strcat(String_Ouput, File_Name);

	FILE *Foutput;
	if ((Foutput = fopen(String_Ouput, "wb")) != NULL)
	{
		fwrite(save_slice, sizeof(unsigned short), sizeX*sizeY, Foutput);
	}
	fclose(Foutput);

	printf("Data Save: %s \n", File_Name);

	delete[] save_slice;                       // Clearing out the buffer memories assigned in this function
	delete[] File_Name;
	delete[] String_Ouput;

}

void writeImageData_slice_FloattoFloat(float* processed, int sizeX, int sizeY, int ith) {

	int	  String_Size = 200;
	char *File_Name = new char[String_Size];
	char *String_Ouput = new char[String_Size];
	float* save_slice = new float[sizeY*sizeX];
	const char* path_s = path.c_str();
	const char* path_Save_s = path_Save.c_str();

	// Writing Voxel Slices
	for (int iny = 0; iny < sizeY; iny++) {
		for (int inx = 0; inx < sizeX; inx++) {
			save_slice[iny *sizeX + inx] = (float)((processed[iny *sizeX + inx]));
		}
	}


	sprintf(File_Name, "/%04d.raw", ith);
	strcpy(String_Ouput, path_Save_s);
	strcat(String_Ouput, "/Result");
	_mkdir(String_Ouput);
	strcat(String_Ouput, File_Name);

	FILE *Foutput;
	if ((Foutput = fopen(String_Ouput, "wb")) != NULL)
	{
		fwrite(save_slice, sizeof(float), sizeX*sizeY, Foutput);
	}
	fclose(Foutput);

	printf("Data Save: %s \n", File_Name);

	delete[] save_slice;                       // Clearing out the buffer memories assigned in this function
	delete[] File_Name;
	delete[] String_Ouput;
}

unsigned short* readImageData_slice(int sizeX, int sizeY, int ith) {


	int	  String_Size = 200;
	char  *String_Input = new char[String_Size];
	char  *File_Name = new char[String_Size];
	const char* path_s = path.c_str();
	//
	unsigned short *load_image = new unsigned short[sizeX * sizeY];
	memset(load_image, 0, sizeof(short)*sizeX * sizeY);
	unsigned short *load_image_Buffer = new unsigned short[sizeX*sizeY];
	memset(load_image_Buffer, 0, sizeof(unsigned short)*sizeX*sizeY);

	sprintf(File_Name, "/CASE %03d_PRE.raw", ith);
	strcpy(String_Input, path_s);
	strcat(String_Input, File_Name);
	FILE *Finput;
	if ((Finput = fopen(String_Input, "rb")) != NULL) {
		fread(load_image_Buffer, sizeof(unsigned short), sizeX*sizeY, Finput);
		for (int iny = 0; iny < sizeY; iny++) {
			for (int inx = 0; inx < sizeX; inx++) {
				load_image[iny *sizeX + inx] = load_image_Buffer[iny *sizeX + inx];
			}
		}
		fclose(Finput);
	}
	else {
		printf("------------- Load Data Error ------------- \n");
		return 0;
	}
	printf("Data load: %s \n", File_Name);
	fclose(Finput);



	delete[] load_image_Buffer;
	delete[] File_Name;
	delete[] String_Input;

	return load_image;
}

std::tuple<unsigned short*, int*> crop_image(unsigned short* img, int * dimensions)
{
	unsigned short max = max_ws(img, dimensions[0] * dimensions[1]);

	unsigned short *img_map = new unsigned short[dimensions[0] * dimensions[1]];
	memset(img_map, 0, sizeof(unsigned short)*dimensions[0] * dimensions[1]);

	for (int i = 0; i < dimensions[0] * dimensions[1]; i++) {
		if (img[i]<max*0.7) {
			img_map[i] = 0;
		}
		else {
			img_map[i] = 1;
		}
	}
	int sizeX = (int)(dimensions[0]);
	int sizeY = (int)(dimensions[1]);
	//find X
	int row_number[2];
	////
	row_number[0] = sizeX;
	row_number[1] = 0;
	////
	short* row_sum = new short[sizeX];
	memset(row_sum, 0, sizeof(short)*sizeX);
	for (int iny = 0; iny < sizeX; iny++) { //start: iny*dimension[1]	//end: iny*dimension[1] + dimension[0]
		for (int inx = iny*sizeY; inx < iny*sizeY + sizeX; inx++) {
			row_sum[iny] = row_sum[iny] + img_map[inx];
			if (row_sum[iny] != 0 && row_sum[iny - 1] == 0 && row_number[0]>iny) {
				row_number[0] = iny;
			}
			else if (row_sum[iny] != 0 && row_sum[iny + 1] == 0 && row_number[1]<iny) {
				row_number[1] = iny;
			}
		}
	}
	//find Y
	int col_number[2];
	////
	col_number[0] = sizeY;
	col_number[1] = 0;
	////
	short* col_sum = new short[sizeY];
	memset(col_sum, 0, sizeof(short)*sizeY);
	for (int inx = 0; inx < sizeY; inx++) {
		for (int iny = 0; iny < sizeX; iny++) {
			col_sum[inx] = col_sum[inx] + img_map[iny*sizeY + inx];
			if (col_sum[inx] != 0 && col_sum[inx - 1] == 0 && col_number[0]>inx) {
				col_number[0] = inx;
			}
			else if (col_sum[inx] != 0 && col_sum[inx + 1] == 0 && col_number[1]<inx) {

				col_number[1] = inx;
			}
		}
	}

	printf("\nrow_min:%d", row_number[0]);
	printf("\nrow_min:%d", row_number[1]);

	printf("\ncol_min:%d", col_number[0]);
	printf("\ncol_min:%d", col_number[1]);

	int* crop_index = new int[1];
	crop_index[0] = col_number[1] - col_number[0] + 1; ////widht
	crop_index[1] = row_number[1] - row_number[0] + 1; ////heigth

	printf("\nwidht:%d", crop_index[0]);
	printf("\nheigth:%d\n", crop_index[1]);

	unsigned short* cropped_image = new unsigned short[crop_index[0] * crop_index[1]];
	for (int iny = 0; iny < crop_index[1]; iny++) {
		for (int inx = 0; inx < crop_index[0]; inx++) {
			cropped_image[iny * crop_index[0] + inx] = (unsigned short)((img[dimensions[1] * (row_number[0] + iny) + col_number[0] + inx]));
		}
	}

	delete[] img_map, img, row_number, row_sum, col_number, col_sum;
	return make_tuple(cropped_image, crop_index);

}

Values Sensor_value(unsigned short* img, int sizeX, int sizeY, int filter_size) {

	unsigned int Max_Intensity, Sensor_V;
	Max_Intensity = img[0];
	for (int i = 0; i < sizeX*sizeY; i++) {
		if (img[i] > Max_Intensity) {
			Max_Intensity = img[i];
		}
	}

	if (Max_Intensity >= unsigned int(40000) && Max_Intensity <= unsigned int(50000)) {
		Sensor_V = unsigned int(60000);
	}
	else if (Max_Intensity >= unsigned int(30000) && Max_Intensity <= unsigned int(40000)) {
		Sensor_V = unsigned short(50000);
	}
	else if (Max_Intensity >= unsigned int(20000) && Max_Intensity <= unsigned int(30000)) {
		Sensor_V = unsigned short(40000);
	}
	else if (Max_Intensity >= unsigned int(0) && Max_Intensity <= unsigned int(20000)) {
		Sensor_V = unsigned int(30000);
	}
	Values Values;

	Values.max = Max_Intensity;
	Values.sensor = Sensor_V;
	Values.padded_size = filter_size / 2;;
	Values.sizeX_pad = sizeX + Values.padded_size * 2;
	Values.sizeY_pad = sizeY + Values.padded_size * 2;

	return Values;
}

unsigned short* DeadpixelCal(unsigned short* img, int sizeX, int sizeY) {

	unsigned short max = max_ws(img, sizeX * sizeY);

	unsigned int *img_map = new unsigned int[sizeX * sizeY];
	memset(img_map, 0, sizeof(unsigned int)*sizeX*sizeY);

	unsigned int count = 0;
	for (int i = 0; i < sizeX*sizeY; i++) {
		if (img[i] > max * 0.991) {
			img_map[count] = i;
			count = count + 1;
		}
	}
	unsigned int *img_max_map = new unsigned int[count - 1];
	memset(img_max_map, 0, sizeof(unsigned int)*count - 1);

	for (int i = 0; i < count; i++) {
		img_max_map[i] = img_map[i];
	}

	unsigned int right_index, left_index, botten_index, top_index;

	for (int i = 0; i < count; i++) {
		if ((img_max_map[i] / sizeX) == sizeX - 1) {//�� �ϴ�
			right_index = sizeX*(sizeX - 2) + mod_fuction(img_max_map[i], sizeX);
		}
		else {
			right_index = sizeX*((img_max_map[i] / sizeX) + 1) + mod_fuction(img_max_map[i], sizeX);
		}
		if ((img_max_map[i] / sizeX) == 0) {// �� ���
			left_index = sizeX + mod_fuction(img_max_map[i], sizeX);
		}
		else {
			left_index = sizeX* ((img_max_map[i] / sizeX) - 1) + mod_fuction(img_max_map[i], sizeX);
		}
		if (mod_fuction(img_max_map[i], sizeY) == 0) {//�� ����
			botten_index = img_max_map[i] - 1;
		}
		else {
			botten_index = img_max_map[i] + 1;
		}
		if (mod_fuction(img_max_map[i], sizeY) == 1) {// �� ����
			top_index = img_max_map[i] + 1;
		}
		else {
			top_index = img_max_map[i] - 1;
		}

		//printf("Dead Pixel:%d\n", img[img_max_map[i]]);
		img[img_max_map[i]] = (img[left_index] + img[right_index] + img[botten_index] + img[top_index])*0.25;
		//printf("Caled Pixel:%d\n", img[img_max_map[i]]);
		//printf("%dth pixel cal complete!\n", img_max_map[i]);

	}
	printf("pixel cal complete!\n");
	delete[] img_map, img_max_map;
	return img;
}

float* Padded_slice_replicate(float* img, int sizeX, int sizeY, int filter_size, Values Values) {

	float *padded_img = new float[Values.sizeX_pad*Values.sizeY_pad];
	memset(padded_img, 0, Values.sizeX_pad*Values.sizeY_pad * sizeof(float));

	for (int j = 0; j < sizeY; j++) {
		for (int i = 0; i < sizeX; i++) {
			padded_img[(j + ((Values.sizeY_pad - sizeY) / 2))*Values.sizeX_pad + i + ((Values.sizeX_pad - sizeX) / 2)] = img[j*sizeX + i];
		}
	}

	for (int j = 0; j < sizeY; j++) {
		for (int i = 0; i < (Values.sizeX_pad - sizeX) / 2; i++) {
			// Left padding
			padded_img[(j + ((Values.sizeY_pad - sizeY) / 2))*Values.sizeX_pad + i] = img[j*sizeX + 0];
			// Right padding
			padded_img[(j + ((Values.sizeY_pad - sizeY) / 2))*Values.sizeX_pad + i + ((Values.sizeX_pad - sizeX) / 2 + sizeX)] = img[j*sizeX + sizeX - 1];
		}
	}
	for (int j = 0; j < (Values.sizeY_pad - sizeY) / 2; j++) {
		for (int i = 0; i < Values.sizeX_pad; i++) {
			// Up padding
			padded_img[j*Values.sizeX_pad + i] = padded_img[(((Values.sizeY_pad - sizeY) / 2) + 0)*Values.sizeX_pad + i];
			// Down padding
			padded_img[(j + ((Values.sizeY_pad - sizeY) / 2) + sizeY)*Values.sizeX_pad + i] = padded_img[(((Values.sizeY_pad - sizeY) / 2) + sizeY - 1)*Values.sizeX_pad + i];
		}
	}

	return padded_img;
}

float* Box_filter(float* padded_img, int sizeX, int sizeY, int filter_size, Values Values) {

	float *Box_result = new float[sizeX * sizeY];
	memset(Box_result, 0, sizeof(float)*sizeX * sizeY);

	float *serch_patch = new float[filter_size * filter_size];
	memset(serch_patch, 0, sizeof(float)*filter_size * filter_size);

	for (int ii = 0; ii < sizeX*sizeY; ii++) {
		for (int inyy = 0; inyy < filter_size; inyy++) {
			for (int inxx = 0; inxx < filter_size; inxx++) {
				serch_patch[inyy *filter_size + inxx] = padded_img[(sizeX + 2 * Values.padded_size)* (Values.padded_size + ii / sizeX + (inyy - Values.padded_size)) + Values.padded_size + ii - (sizeX)*(ii / sizeX) + (inxx - Values.padded_size)];
			}
		}

		float means = mean(serch_patch, filter_size * filter_size);
		Box_result[ii] = means;
		//float eps = std::numeric_limits<float>::epsilon();
		//Box_result[ii] = (padded_img[(sizeX + 2 * Values.padded_size)* (Values.padded_size + ii / sizeX) + Values.padded_size + ii - (sizeX)*(ii / sizeX)] - Tx + eps) / (1 - Tx + eps);
	}

	delete[] serch_patch;
	return Box_result;

}


int main()
{
	int Dimension[1], ith, filter_size;
	float gamma, alpha, epsilon;
	string crop_option;
	Dimension[0] = 3072;	////input: x size ũ��
	Dimension[1] = 3072;	////input: y size ũ��
	ith = 12;					////input: �� ��° ������ ����?
	filter_size = 3;		////input: local enhancement�� ���� ������
	crop_option = "on";		////input: activation ���� crop ����

	epsilon = 0.0002;

							////Data loading �κ�
	unsigned short *O_Data = readImageData_slice(Dimension[0], Dimension[1], ith);
	////Dead pixel calibration �κ�
	DeadpixelCal(O_Data, Dimension[0], Dimension[1]);

	unsigned short * Input_Data_crop;
	int * crop_index;
	////activation ���� crop �κ�
	if (crop_option == "on") {
		tie(Input_Data_crop, crop_index) = crop_image(O_Data, Dimension);
		writeImageData_slice(Input_Data_crop, crop_index[0], crop_index[1], string("Cropped"), ith);
		Dimension[0] = crop_index[0];
		Dimension[1] = crop_index[1];
	}
	else {
		Input_Data_crop = O_Data;
	}

	unsigned short *Input_Data = new unsigned short[Dimension[0] * Dimension[1]];
	memset(Input_Data, 0, sizeof(unsigned short)*(Dimension[0] * Dimension[1]));
	for (int i = 0; i < Dimension[0] * Dimension[1]; i++) {
		Input_Data[i] = Input_Data_crop[i];
	}
	////����ϴ� parameter ����
	Values Values;
	Values = Sensor_value(Input_Data, Dimension[0], Dimension[1], filter_size);

	float *img_reverse = new float[Dimension[0] * Dimension[1]];
	memset(img_reverse, 0, sizeof(float)*Dimension[0] * Dimension[1]);
	////���� ���� ��Ʈ
	for (int i = 0; i < Dimension[0] * Dimension[1]; i++) {
		img_reverse[i] = float(Input_Data[i]) / Values.sensor;
		img_reverse[i] = 1 - img_reverse[i];
	}



	//����ϴ� map ����
	float *padded_img_P = Padded_slice_replicate(img_reverse, Dimension[0], Dimension[1], filter_size, Values);
	float *padded_img_I = Padded_slice_replicate(img_reverse, Dimension[0], Dimension[1], filter_size, Values);

	float *Mu_P = Box_filter(padded_img_P, Dimension[0], Dimension[1], filter_size, Values);
	float *Mu_I = Box_filter(padded_img_I, Dimension[0], Dimension[1], filter_size, Values);

	float *padded_img_I_P = new float[Values.sizeX_pad*Values.sizeY_pad];
	memset(padded_img_I_P, 0, sizeof(float)*(Values.sizeX_pad*Values.sizeY_pad));
	for (int i = 0; i < Values.sizeX_pad*Values.sizeY_pad; i++) {
		padded_img_I_P[i] = padded_img_I[i] *padded_img_P[i];
	}
	float *Cross_I = Box_filter(padded_img_I_P, Dimension[0], Dimension[1], filter_size, Values);

	float *padded_img_I_I = new float[Values.sizeX_pad*Values.sizeY_pad];
	memset(padded_img_I_I, 0, sizeof(float)*(Values.sizeX_pad*Values.sizeY_pad));
	for (int i = 0; i < Values.sizeX_pad*Values.sizeY_pad; i++) {
		padded_img_I_I[i] = padded_img_I[i] * padded_img_I[i];
	}
	float *Sig_I = Box_filter(padded_img_I_I, Dimension[0], Dimension[1], filter_size, Values);

	//*padded_img_I_P ** padded_img_I_P �̰� �Ǹ� for�� �Ⱦ��ٵ�...

	float *Var_I = new float[Dimension[0]* Dimension[1]];
	memset(Var_I, 0, sizeof(float)*(Dimension[0]* Dimension[1]));
	for (int i = 0; i < Dimension[0]* Dimension[1]; i++) {
		Var_I[i] = Sig_I[i] - pow(Mu_I[i],2);
	}

	///guided filter���� �����ϴ� map ��� �ٸ� ���� ����ϰ� �ʹٸ�,
	///aaa�� edge map���� ��ü �ϸ� �˴ϴ�
	float *aaa = new float[Dimension[0]* Dimension[1]];
	float *bbb = new float[Dimension[0]* Dimension[1]];

	for (int i = 0; i < Dimension[0]* Dimension[1]; i++) {
		aaa[i] = (Cross_I[i] - (Mu_I[i]* Mu_P[i])) / (epsilon + Var_I[i]);
	}

	for (int i = 0; i < Dimension[0]* Dimension[1]; i++) {
		bbb[i] = Mu_P[i] - aaa[i]* Mu_I[i];
	}

	float *padded_aaa = Padded_slice_replicate(aaa, Dimension[0], Dimension[1], filter_size, Values);
	float *padded_bbb = Padded_slice_replicate(bbb, Dimension[0], Dimension[1], filter_size, Values);


	float *Mu_A = Box_filter(padded_aaa, Dimension[0], Dimension[1], filter_size, Values);
	float *Mu_B = Box_filter(padded_bbb, Dimension[0], Dimension[1], filter_size, Values);


	float *Result = new float[Dimension[0]* Dimension[1]];
	memset(Result, 0, sizeof(float)*(Dimension[0]* Dimension[1]));
	for (int i = 0; i < Dimension[0]* Dimension[1]; i++) {
		Result[i] = Mu_A[i]* img_reverse[i] + Mu_B[i];
	}

	////14bit saving
	for (int i = 0; i < Dimension[0] * Dimension[1]; i++) {
		img_reverse[i] = img_reverse[i] * pow(2, 14);
	}

	for (int i = 0; i < Dimension[0] * Dimension[1]; i++) {
		Result[i] = Result[i] * pow(2, 14);
	}


	//writeImageData_slice(Input_Data, Dimension[0], Dimension[1], string("Ori_deadcal"), ith);
	writeImageData_slice_FloattoShort(img_reverse, Dimension[0], Dimension[1], string("Input"), ith);
	writeImageData_slice_FloattoShort(Result, Dimension[0], Dimension[1], string("Result"), ith);
	printf("end");
}



