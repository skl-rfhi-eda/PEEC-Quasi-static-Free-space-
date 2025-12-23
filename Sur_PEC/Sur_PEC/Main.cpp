#include "Mod_Read_File.h"
#include "Mod_Connection.h"
#include "Mod_Element_Built.h"
#include "Mod_Solver_Freq.h"
#include "Mod_Solver_Time.h"


int main(){

	Ini_Gaussion();
	Read_File();
	Creat_Connection();
	Element_Built();

	if (Solver_SET == 0)
		Freq_Solver();
	else if (Solver_SET == 1)
		Time_Solver();
	system("pause");
}