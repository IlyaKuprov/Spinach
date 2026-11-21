#include "../StdAfx.h"
#include "isotopes.h"

/**
*\file isotopes.cpp
*/

/**
*  @brief Main function for creating static isotopes in the Isotopes collection in order
*/
void Isotopes::Init_Isotopes()
{
	Isotopes_Collection= gcnew ArrayList(); // initialization of array
	float rev=1.0f/255.0f;                  //Reverse of 255, for colours

	//Create each isotope with the respected properties and then add to the collection in order
	Isotope^ isotope;
	//10
	isotope=gcnew Isotope("e",   0,   0,  0,   0.1f, 0.1f, 0.1f, 0.5f, 0.0f, -1.760859770e11);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("H",   1,   1,  53,  204.0f*rev, 204.0f*rev, 204.0f*rev, 0.5f, 0.0f, 26.7522205e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("H",   1,   2,  53,  204.0f*rev, 204.0f*rev, 204.0f*rev, 1.0f, 0.00286f, 4.10662914e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("H",   1,   3,  53,  204.0f*rev, 204.0f*rev, 204.0f*rev, 0.5f, 0.0f, 28.5350e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("He",  2,   3,  31,  216.0f*rev, 1.0f, 1.0f, 0.5f, 0.0f, -20.3801e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Li",  3,   6,  145, 204.0f*rev, 124.0f*rev, 1.0f, 1.0f, -0.00083f, 3.9371e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Li",  3,   7,  145, 204.0f*rev, 124.0f*rev, 1.0f, 1.5f, -0.04f, 10.3976e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Be",  4,   9,  105, 204.0f*rev, 1.0f, 0.0f, 1.5f, 0.053f, -3.7606e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("B",   5,   10, 85,  1.0f, 181.0f*rev, 181.0f*rev, 3.0f, 0.0847f, 2.8746786e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("B",   5,   11, 85,  1.0f, 181.0f*rev, 181.0f*rev, 1.5f, 0.0407f, 8.5847e7);   Isotopes_Collection->Add(isotope);

	//20
	isotope=gcnew Isotope("C",   6,   13,  70,  102.0f*rev, 102.0f*rev, 102.0f*rev, 0.5f, 0.0f, 6.728284e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("N",   7,   14,  65,  25.0f*rev, 25.0f*rev, 229.0f*rev, 1.0f, 0.02f, 1.9337792e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("N",   7,   15,  65,  25.0f*rev, 25.0f*rev, 229.0f*rev, 0.5f, 0.0f, -2.71261804e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("O",   8,   17,  60,  229.0f*rev, 0.0f, 0.0f, 2.5f, -0.02578f, -3.62808e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("F",   9,   19,  50,  178.0f*rev, 1.0f, 1.0f, 0.5f, 0.0f, 2.518148e8);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ne",  10,  21,  38,  175.0f*rev, 226.0f*rev, 244.0f*rev, 1.5f, 0.103f, -2.11308e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Na",  11,  23,  180, 170.0f*rev, 91.0f*rev, 242.0f*rev, 1.5f, 0.101f, 7.0704e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Mg",  12,  25,  150, 178.0f*rev, 204.0f*rev, 0.0f, 2.5f, 1.25f, -1.6389e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Al",  13,  27,  125, 209.0f*rev, 165.0f*rev, 165.0f*rev, 2.5f, 0.1402f, 6.971e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Si",  14,  29,  110, 127.0f*rev, 153.0f*rev, 153.0f*rev, 0.5f, 0.0f, -5.314e7);   Isotopes_Collection->Add(isotope);

	//30
	isotope=gcnew Isotope("P",   15,  31,  100, 1.0f, 127.0f*rev, 0.0f, 0.5f, 0.0f, 10.8394e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("S",   16,  33,  100, 1.0f, 198.0f*rev, 40.0f*rev, 1.5f, -0.064f, 2.0557e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cl",  17,  35,  100, 25.0f*rev, 239.0f*rev, 25.0f*rev, 1.5f, 0.0f, 2.6242e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cl",  17,  37,  100, 25.0f*rev, 239.0f*rev, 25.0f*rev, 1.5f, -0.06493f, 2.1844e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ar", 18, 40, 71, 0.5f, 0.8f, 0.8f, 0.0f, 0.0f, 0.0); Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("K",   19,  39,  220, 142.0f*rev, 63.0f*rev, 211.0f*rev, 1.5f, 0.049f, 1.2499e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("K",   19,  41,  220, 142.0f*rev, 63.0f*rev, 211.0f*rev, 1.5f, 0.06f, 0.6861e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ca",  20,  43,  180, 153.0f*rev, 153.0f*rev, 0.0f, 3.5f, -0.043f, -1.8028e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sc",  21,  45,  160, 229.0f*rev, 229.0f*rev, 226.0f*rev, 3.5f, -0.22f, 6.5088e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ti",  22,  47,  140, 191.0f*rev, 193.0f*rev, 198.0f*rev, 2.5f, 0.3f, -1.5106e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ti",  22,  49,  140, 191.0f*rev, 193.0f*rev, 198.0f*rev, 3.5f, 0.24f,-1.5110e7);   Isotopes_Collection->Add(isotope);

	//40
	isotope=gcnew Isotope("V",   23,  50, 135,  165.0f*rev, 165.0f*rev, 170.0f*rev, 6.0f, 0.21f, 2.6721e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("V",   23,  51, 135,  165.0f*rev, 165.0f*rev, 170.0f*rev, 3.5f, -0.043f, 7.0492e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cr",  24,  53, 166,  137.0f*rev, 153.0f*rev, 198.0f*rev, 1.5f, -0.15f, -1.5077e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Mn",  25,  55, 140,  155.0f*rev, 122.0f*rev, 198.0f*rev, 2.5f, 0.33f, 6.6453e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Fe",  26,  57, 140,  127.0f*rev, 122.0f*rev, 198.0f*rev, 0.5f, 0.0f, 0.8687e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Co",  27,  59, 135,  91.0f*rev, 109.0f*rev, 1.0f, 3.5f, 0.35f, 6.3015e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ni",  28,  61, 135,  91.0f*rev, 122.0f*rev, 193.0f*rev, 1.5f, 0.0162f, -2.3943e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cu",  29,  63, 135,  1.0f, 122.0f*rev, 96.0f*rev, 1.5f, -0.211f, 7.1088e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cu",  29,  65, 135,  1.0f, 122.0f*rev, 96.0f*rev, 1.5f, -0.195f, 7.6104e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Zn",  30,  67, 135,  124.0f*rev, 127.0f*rev, 175.0f*rev, 2.5f, 0.15f, 1.6778e7);   Isotopes_Collection->Add(isotope);

	//50
	isotope=gcnew Isotope("Ga",  31,  69, 130,  193.0f*rev, 142.0f*rev, 142.0f*rev, 1.5f, 2.01659f, 6.4389e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ga",  31,  71, 130,  193.0f*rev, 142.0f*rev, 142.0f*rev, 1.5f, 0.106f, 8.1812e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ge",  32,  73, 125,  102.0f*rev, 142.0f*rev, 142.0f*rev, 4.5f, -0.17f, -0.93660e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("As",  33,  75, 115,  188.0f*rev, 127.0f*rev, 226.0f*rev, 1.5f, 0.3f, 4.5961e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Se",  34,  77, 115,  1.0f, 160.0f*rev, 0.0f, 0.5f, 0.0f, 5.1214e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Br",  35,  79, 115,  165.0f*rev, 33.0f*rev, 33.0f*rev, 1.5f, 0.305f, 6.7256e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Br",  35,  81, 115,  165.0f*rev, 33.0f*rev, 33.0f*rev, 4.5f, 0.254f, 7.249776e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Kr",  36,  83,  88,  91.0f*rev, 186.0f*rev, 209.0f*rev, 3.5f, 0.26f, -1.0331e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Rb",  37,  85, 235,  112.0f*rev, 45.0f*rev, 175.0f*rev, 2.5f, 0.23f, 2.5923e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Rb",  37,  87, 235,  112.0f*rev, 45.0f*rev, 175.0f*rev, 1.5f, 0.132f, 8.7851e7);   Isotopes_Collection->Add(isotope);

	//60
	isotope=gcnew Isotope("Sr",  38,  87,  200, 127.0f*rev, 102.0f*rev, 0.0f, 4.5f, 0.34f, -1.1635e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Y",   39,  89,  180, 147.0f*rev, 252.0f*rev, 1.0f, 0.5f, 0.0f, -1.3163e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Zr",  40,  91,  155, 147.0f*rev, 224.0f*rev, 224.0f*rev, 2.5f, -0.206f, -2.4975e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Nb",  41,  93,  145, 114.0f*rev, 193.0f*rev, 201.0f*rev, 4.5f, -0.32f, 6.5674e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Mo",  42,  95,  145, 84.0f*rev, 181.0f*rev, 181.0f*rev, 2.5f, -0.022f, -1.7514e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Mo",  42,  97,  145, 84.0f*rev, 181.0f*rev, 181.0f*rev, 2.5f, 0.255f, -1.7884e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Tc",  43,  99,  135, 58.0f*rev, 158.0f*rev, 168.0f*rev, 4.5f, -0.129f, 6.0503e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ru",  44,  99,  130, 35.0f*rev, 142.0f*rev, 150.0f*rev, 2.5f, 0.079f, -1.2286e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ru",  44,  101, 130, 35.0f*rev, 142.0f*rev, 150.0f*rev, 2.5f, 0.46f, -1.3773e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Rh",  45,  103, 135, 10.0f*rev, 124.0f*rev, 140.0f*rev, 0.5f, 0.0f, -0.8468e7);   Isotopes_Collection->Add(isotope);

	//70
	isotope=gcnew Isotope("Pd",  46,  105, 140, 0.0f, 104.0f*rev, 132.0f*rev, 2.5f, 0.660f, -1.2305e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ag",  47,  107, 160, 153.0f*rev, 198.0f*rev, 1.0f, 0.5f, 0.0f, -1.0878e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ag",  47,  109, 160, 153.0f*rev, 198.0f*rev, 1.0f, 0.5f, 0.0f, -1.2519e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cd",  48,  111, 155, 1.0f, 216.0f*rev, 142.0f*rev, 0.5f, 0.0f, -5.7046e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cd",  48,  113, 155, 1.0f, 216.0f*rev, 142.0f*rev, 0.5f, 0.0f, -5.9609e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("In",  49,  113, 155, 165.0f*rev, 117.0f*rev, 114.0f*rev, 4.5f, 0.08f, 5.8845e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("In",  49,  115, 155, 165.0f*rev, 117.0f*rev, 114.0f*rev, 4.5f, 0.81f, 5.8971e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sn",  50,  117, 145, 102.0f*rev, 127.0f*rev, 127.0f*rev, 0.5f, 0.0f, -9.589e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sn",  50,  119, 145, 102.0f*rev, 127.0f*rev, 127.0f*rev, 0.5f, 0.0f, -10.0318e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sb",  51,  121, 145, 158.0f*rev, 99.0f*rev, 181.0f*rev, 2.5f, -0.36f, 6.4442e7);   Isotopes_Collection->Add(isotope);

	//80
	isotope=gcnew Isotope("Sb",  51,  123, 145, 158.0f*rev, 99.0f*rev, 181.0f*rev, 3.5f, -0.49f, 3.4904e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Te",  52,  123, 140, 211.0f*rev, 122.0f*rev, 0.0f, 0.5f, 0.0f, -7.0576e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Te",  52,  125, 140, 211.0f*rev, 122.0f*rev, 0.0f, 0.5f, 0.0f, -8.5087e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("I",   53,  127, 140, 147.0f*rev, 0.0f, 147.0f*rev, 2.5f, 0.689f, 5.3896e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Xe",  54,  129, 108, 66.0f*rev, 158.0f*rev, 175.0f*rev, 0.5f, 0.0f, -7.4521e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Xe",  54,  131, 108, 66.0f*rev, 158.0f*rev, 175.0f*rev, 1.5f, -0.116f, 2.2091e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Cs",  55,  133, 260, 86.0f*rev, 22.0f*rev, 142.0f*rev, 3.5f, -0.00371f, 3.5339e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ba",  56,  135, 215, 102.0f*rev, 51.0f*rev, 0.0f, 1.5f, 0.16f, 2.657e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ba",  56,  137, 215, 102.0f*rev, 51.0f*rev, 0.0f, 1.5f, 0.245f, 2.973e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("La",  57,  139, 195, 112.0f*rev, 221.0f*rev, 1.0f, 3.5f, 0.2f, 3.8085e7);   Isotopes_Collection->Add(isotope);

	//90
	isotope=gcnew Isotope("Pr",  59,  141, 185, 216.0f*rev, 1.0f, 198.0f*rev, 2.5f, -0.077f, 7.765e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Nd",  60,  143, 185, 198.0f*rev, 1.0f, 198.0f*rev, 3.5f, -0.61f, -1.474e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Nd",  60,  145, 185, 198.0f*rev, 1.0f, 198.0f*rev, 3.5f, -0.314f, -0.913e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sm",  62,  147, 185, 142.0f*rev, 1.0f, 198.0f*rev, 3.5f, -0.27f,-1.1124e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sm",  62,  149, 185, 142.0f*rev, 1.0f, 198.0f*rev, 3.5f, -0.078f, -0.9175e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Eu",  63,  151, 185, 96.0f*rev, 1.0f, 198.0f*rev, 2.5f, 0.83f,  6.5477e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Eu",  63,  153, 185, 96.0f*rev, 1.0f, 198.0f*rev, 2.5f, 2.22f, 2.9371e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Gd",  64,  155, 180, 68.0f*rev, 1.0f, 198.0f*rev, 1.5f, 1.27f, -0.8273e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Gd",  64,  157, 180, 68.0f*rev, 1.0f, 198.0f*rev, 1.5f, 1.36f, -1.0792e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Tb",  65,  159, 175, 48.0f*rev, 1.0f, 198.0f*rev, 1.5f, 1.432f, 6.4306e7);   Isotopes_Collection->Add(isotope);

	//100
	isotope=gcnew Isotope("Dy",  66,  161, 175, 30.0f*rev, 1.0f, 181.0f*rev, 2.5f, 2.51f, -0.9206e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Dy",  66,  163, 175, 30.0f*rev, 1.0f, 181.0f*rev, 2.5f, 2.65f, 1.275e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ho",  67,  165, 175, 0.0f, 1.0f, 181.0f*rev, 2.5f, 3.58f, 5.487e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Er",  68,  167, 175, 0.0f, 229.0f*rev, 117.0f*rev, 3.5f, 3.57f, -0.7752e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Tm",  69,  169, 175, 0.0f, 211.0f*rev, 81.0f*rev, 0.5f, 0.0f,  -2.1376e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Yb",  70,  171, 175, 0.0f, 191.0f*rev, 56.0f*rev, 0.5f, 0.0f, 4.7348e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Yb",  70,  173, 175, 0.0f, 191.0f*rev, 56.0f*rev, 2.5f, 2.8f, -1.3025e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Lu",  71,  175, 175, 0.0f, 170.0f*rev, 35.0f*rev, 3.5f, 3.49f, 3.0589e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Lu",  71,  176, 175, 0.0f, 170.0f*rev, 35.0f*rev, 3.5f, 4.92f, 2.1684e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Hf",  72,  177, 155, 76.0f*rev, 193.0f*rev, 1.0f, 3.5f, 3.37f, 1.086e7);   Isotopes_Collection->Add(isotope);

	//110
	isotope=gcnew Isotope("Hf",  72,  179, 155, 76.0f*rev, 193.0f*rev, 1.0f, 4.5f, 3.79f, -0.6821e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ta",  73,  181, 145, 76.0f*rev, 165.0f*rev, 1.0f, 3.5f, 3.17f, 3.2445e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("W",   74,  183, 135, 38.0f*rev, 147.0f*rev, 214.0f*rev, 0.5f, 0.0f, 1.1283e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Re",  75,  185, 135, 38.0f*rev, 124.0f*rev, 170.0f*rev, 2.5f, 2.18f, 6.1057e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Re",  75,  187, 135, 38.0f*rev, 124.0f*rev, 170.0f*rev, 2.5f, 2.07f, 6.1682e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Os",  76,  187, 130, 38.0f*rev, 102.0f*rev, 150.0f*rev, 0.5f, 0.0f, 0.6193e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Os",  76,  189, 130, 38.0f*rev, 102.0f*rev, 150.0f*rev, 1.5f, 0.86f, 2.1072e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ir",  77,  191, 135, 22.0f*rev, 84.0f*rev, 135.0f*rev, 1.5f, 2.5f, 0.4665e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Ir",  77,  193, 135, 22.0f*rev, 84.0f*rev, 135.0f*rev, 1.5f, 0.47f, 0.508e7);   Isotopes_Collection->Add(isotope);

	//120
	isotope=gcnew Isotope("Pt",  78,  195, 135, 22.0f*rev, 91.0f*rev, 142.0f*rev, 0.5f, 0.0f, 5.8383e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Au",  79,  197, 135, 1.0f, 209.0f*rev, 35.0f*rev, 1.5f, 0.547f, 0.4692e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Hg",  80,  199, 150, 181.0f*rev, 181.0f*rev, 193.0f*rev, 0.5f, 0.0f, 4.8458e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Hg",  80,  201, 150, 181.0f*rev, 181.0f*rev, 193.0f*rev, 1.5f, 0.38f, -1.7888e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Tl",  81,  203, 190, 165.0f*rev, 84.0f*rev, 76.0f*rev, 0.5f, 0.0f, 15.5394e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Tl",  81,  205, 190, 165.0f*rev, 84.0f*rev, 76.0f*rev, 0.5f, 0.0f, 15.6922e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Pb",  82,  207, 180, 86.0f*rev, 89.0f*rev, 96.0f*rev, 0.5f, 0.0f, 5.6264e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Bi",  83,  209, 160, 158.0f*rev, 79.0f*rev, 181.0f*rev, 4.5f, -0.37f, 4.3752e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("U",   92,  235, 175, 0.0f, 142.0f*rev, 1.0f, 3.5f, 4.936f, -0.4926e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Pu",  94,  239, 175, 0.0f, 107.0f*rev, 242.0f*rev, 0.5f, 0.0f, 0.972e7);   Isotopes_Collection->Add(isotope);


	isotope=gcnew Isotope("Am",  95,  243, 175, 84.0f*rev, 91.0f*rev, 242.0f*rev, 2.5f, 2.86f, 1.54e7);   Isotopes_Collection->Add(isotope);
	isotope=gcnew Isotope("Sg", 106, 266, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0); Isotopes_Collection->Add(isotope);

	// Extra elements required by Cartesian imports. Choose the most abundant ground-state nuclide,
	// otherwise the longest lived in etc/isotopes.tsv; Mc through Og use IUPAC (4 May 2022).
	// Nuclear constants follow etc/build_isotopes.m. Unconfirmed or unavailable values stay NaN.
	// Radius and grey colour are display defaults only.
	const float unknown=System::Single::NaN;
	Isotopes_Collection->Add(gcnew Isotope("Ce", 58, 140, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0));
	Isotopes_Collection->Add(gcnew Isotope("Pm", 61, 145, 175, 0.5f, 0.5f, 0.5f, 2.5f, 0.21f, 72799131.89));
	Isotopes_Collection->Add(gcnew Isotope("Po", 84, 209, 175, 0.5f, 0.5f, 0.5f, 0.5f, 0.0f, unknown));
	Isotopes_Collection->Add(gcnew Isotope("At", 85, 210, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Rn", 86, 222, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0));
	Isotopes_Collection->Add(gcnew Isotope("Fr", 87, 223, 175, 0.5f, 0.5f, 0.5f, 1.5f, 1.17f, 37357449.26));
	Isotopes_Collection->Add(gcnew Isotope("Ra", 88, 226, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0));
	Isotopes_Collection->Add(gcnew Isotope("Ac", 89, 227, 175, 0.5f, 0.5f, 0.5f, 1.5f, 1.7f, 35122388.19));
	Isotopes_Collection->Add(gcnew Isotope("Th", 90, 232, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0));
	Isotopes_Collection->Add(gcnew Isotope("Pa", 91, 231, 175, 0.5f, 0.5f, 0.5f, 1.5f, -1.72f, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Np", 93, 237, 175, 0.5f, 0.5f, 0.5f, 2.5f, 3.886f, 60155072.14));
	Isotopes_Collection->Add(gcnew Isotope("Cm", 96, 247, 175, 0.5f, 0.5f, 0.5f, 4.5f, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Bk", 97, 247, 175, 0.5f, 0.5f, 0.5f, 1.5f, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Cf", 98, 251, 175, 0.5f, 0.5f, 0.5f, 0.5f, 0.0f, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Es", 99, 252, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Fm", 100, 257, 175, 0.5f, 0.5f, 0.5f, 4.5f, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Md", 101, 258, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("No", 102, 259, 175, 0.5f, 0.5f, 0.5f, 4.5f, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Lr", 103, 266, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Rf", 104, 267, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Db", 105, 268, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Bh", 107, 270, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Hs", 108, 269, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Mt", 109, 277, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Ds", 110, 282, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0));
	Isotopes_Collection->Add(gcnew Isotope("Rg", 111, 282, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Cn", 112, 285, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Nh", 113, 286, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Fl", 114, 290, 175, 0.5f, 0.5f, 0.5f, 0.0f, 0.0f, 0.0));
	Isotopes_Collection->Add(gcnew Isotope("Mc", 115, 290, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Lv", 116, 293, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Ts", 117, 294, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
	Isotopes_Collection->Add(gcnew Isotope("Og", 118, 294, 175, 0.5f, 0.5f, 0.5f, unknown, unknown, unknown));
};

/**
* @brief Function for finding the isotope when one of the element, atom number or mass is given.
* @param element Name of the element.
* @param atomNumber Atomic number of the element.
* @param mass Mass number of the element.
* @result Returns the isotope object found from the above properties.
*/
Isotope^ Isotopes::FindIsotope(System::String^ element, int atomNumber, int mass) 
{
	//Initialize isotopes collection
	Init_Isotopes();

	//Retain the single magnetic-isotope default used by coordinate importers.
	Isotope^ only=nullptr;
	int count=0;
	if(element!=nullptr)
		for each(Isotope^ isotope in Isotopes_Collection)
			if(isotope->Element==element) {only=isotope;count++;}
	if(count==1) return only;

	//Select the first matching isotope when no mass is specified.
	for each(Isotope^ isotope in Isotopes_Collection)
	{
		if(element!=nullptr && isotope->Element!=element) continue;
		if(element==nullptr && atomNumber!=0 && isotope->AtomNumber!=atomNumber) continue;
		if(mass!=0 && isotope->Mass!=mass) continue;
		return isotope;
	}
	throw gcnew System::String("Unknown isotope: "+element+" ("+atomNumber+", "+mass+").");
}

/**
* @brief Function for finding the isotope of the lement but with higher spin.
* @param isotope Isotope object.
* @result Returns the isotope object found with the same element but higher spin.
*/
Isotope^ Isotopes::FindIsotopeWithHigherSpin(Isotope^ isotope) 
{
	//Initialize isotopes collection
	Init_Isotopes();

	//Quadrupolar interactions need spin greater than one half, at any mass.
	for each(Isotope^ candidate in Isotopes_Collection)
		if(candidate->Element==isotope->Element && candidate->Spin>0.5f) return candidate;
	return isotope;
}

/**
* @brief Function for returning all the possible mass numbers of an element. Used for dropmenu.
* @param element Element name.
* @result Returns an array with all the possible masses.
*/
cli::array<System::String^>^ Isotopes::getAllIsotopeMasses(System::String ^element)
{
	Init_Isotopes();
	ArrayList^ masses=gcnew ArrayList();
	for each(Isotope^ isotope in Isotopes_Collection)
		if(isotope->Element==element) masses->Add(isotope->Mass.ToString());
	return safe_cast<cli::array<System::String^>^>(masses->ToArray(System::String::typeid));
}
