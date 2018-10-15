program SampSim_vES;

{$APPTYPE CONSOLE}

{$R *.res}

uses
  System.SysUtils,
  Math,
  Classes,
  LibFiles,
  LibText,
  LibTypes,
  LibjbDBF,
  LibProbDistrib;

Type
  Tsampling = Record
    MatObsLat0: TDblMatrix;
    MatEstLat0: TDblMatrix;
    MatCoverageLat0: TDblMatrix;
    MatGapsLat0: TDblMatrix;
    MatObsCell0: TDblMatrix;
    MatEstCell0: TDblMatrix;
   End;

var
  // Input Matrices
  neighMat, CellLatMat: TDblMatrix;

  // Model Parameters
  Total_ssp, TotalSimu: Integer;

  // Loop control
  Species, nSimu, Scont: Integer; // For model parameters
  idx: Integer; // Generic

  // TLists & pointers
  Cell_ID, SP_ID: TIntVector;
  TRichLat, TObsLat, TRichCell, TObsCell, TSampCell: Array of TList;
  pSP: ^Integer;

  // Generic variables of the equalSamplingApply and CountGapsLat procedures
  // A small modification to avoid resize these arrays at each sampling simulation, which makes the code slow
  sp, vec: TIntVector;
  MatGaps: TDblMatrix;

  // User preference
  Model, inputDir, outputDir: String;
  maxSamp, nSamp: Integer;

  // Output Matrices
  MatRichLat, MatObsLat, MatEstLat, MatCoverageLat, MatGapsLat: TDblMatrix;
  MatRichCell, MatObsCell, MatEstCell: TDblMatrix;
  TEffort: Array of Tsampling;



procedure SpreadingDye;
 var
 LimRange, Founder: Integer;
 i, x: Integer;
 TCell_Ocup, TCell_Neigh: TList;
 pOcup, pNeigh: ^Integer;

 begin

   // Find a randon range size for this species (Note that the species must colonize at least 1 cell)
   //LimRange:= 1 + Random(Length(neighMat)-1);
   LimRange:= Ceil(RndBeta(3, 9)*Length(neighMat));
   // Find a randon cell to start the range filling
   Founder:= Random(Length(neighMat));
   //Founder:= Round(RndBeta(3, 3)*Length(neighMat));

   // Save the ID of the species in the TList of this first cell (Remember to use 'SP_ID' as reference to pointer because 'Species' change through the simulation)
   // This step is repeated to each new cell occupied
   pSP:=@SP_ID[Species];
   TRichCell[Founder].Add(pSP);

   // Create a TList of occupied and neighbors cells
   TCell_Ocup:= TList.Create;
   TCell_Neigh:= TList.Create;
   // Save the starting point in the TList of occupied cells
   pOcup:=@Cell_ID[Founder];
   TCell_Ocup.Add(pOcup);


   // Start the range filling process
   While TCell_Ocup.Count <> LimRange do
   begin

     // First, find the neighbors cells
     // Here I will run all the lines of the neighborhood matrix to find the neighbors of the target cell (in the column)
     // The target cell is pointed by "pOcup", which is the Founder in the begining and randomly sampled after that
     // That is, I will always look for the neighbors of the last cell added to the list
     For i:= 0 to Length(neighMat)-1 do
     begin

       // If the cell is a neighbor, I will test if it can be occupied
       If (neighMat[i,pOcup^] = 1) then
       begin

         // It is important use Cell_ID[] because it is my fixed reference
         // Note that I will use the memory address to test if the cell is already in the TLists
         pNeigh:= @Cell_ID[i];

         // Test if the cell (memory address) is already in the neighbors TList
         If TCell_Neigh.IndexOf(pNeigh) < 0 then
         begin

           // Test if the cell is already occupied
           If TCell_Ocup.IndexOf(pNeigh) < 0 then
           begin

             // If the neighbor cell is not occupied and it is not in the neighbors TList, add to the neighbors TList
             // This way I will have a list only with potential cells for occupation
             TCell_Neigh.Add(pNeigh);
           end;

         end;

       end;

     end;

     // Now sample a neighbor cell to be occuppied
     x:= Random(TCell_Neigh.Count);

     // Again, I have to use the pointer because I'm working with memory address
     // Here I am transfering the address in the position 'x' to the pointer to add it in the list of occuppied cells
     pOcup:= TCell_Neigh[x];
     TCell_Ocup.Add(pOcup);

     // Delete the address of the now occupied cell from the neighbor TList
     TCell_Neigh.Delete(x);

     //Inform the presence of this species in the new occupied cell
     TRichCell[pOcup^].Add(pSP);

   end;

   // The species geographical range was created and the information is saved in the TRichCell TList
   // Finish the other TList of cells
   TCell_Ocup.Free;
   TCell_Neigh.Free;
 end;
procedure ScatteredRange;
 var
 LimRange, Founder: Integer;
 x: Integer;
 TCell_Ocup: TList;
 pOcup: ^Integer;

 begin

   // Find a randon range size for this species (Note that the species must colonize at least 1 cell)
   //LimRange:= 1 + Random(Length(neighMat)-1);
   LimRange:= Round(RndBeta(3, 9)*Length(neighMat));
   // Find a randon cell to start the range filling
   Founder:= Random(Length(neighMat));
   //Founder:= Round(RndBeta(3, 3)*Length(neighMat));

   // Save the ID of the species in the TList of this first cell (Remember to use 'SP_ID' as reference to pointer because 'Species' change through the simulation)
   // This step is repeated to each new cell occupied
   pSP:=@SP_ID[Species];
   TRichCell[Founder].Add(pSP);

   // Create a TList of occupied and neighbors cells
   TCell_Ocup:= TList.Create;
   // Save the starting point in the TList of occupied cells
   pOcup:=@Cell_ID[Founder];
   TCell_Ocup.Add(pOcup);


   // Start the range filling process
   While TCell_Ocup.Count <> LimRange do
   begin

     // Find a randon cell to continue the range filling
     x:= Random(Length(neighMat));
     pOcup:=@Cell_ID[x];

     // Test if the cell is already occupied
     If TCell_Ocup.IndexOf(pOcup) < 0 then
     begin

       // If the cell is not occupied, add to the neighbors TList
       TCell_Ocup.Add(pOcup);

       //Inform the presence of this species in the new occupied cell
       TRichCell[pOcup^].Add(pSP);
     end;

   end;

   // The species geographical range was created and the information is saved in the TRichCell TList
   // Finish the other TList of cells
   TCell_Ocup.Free;
 end;

procedure CountRichLat;
 var
 i, j, cell: Integer;

 begin

   // For each latitude I will calculate the real/known species richness
   For i:=0 to Length(CellLatMat)-1 do
   begin

     // Start a TList for each latitude to save the ID of all species present
     TRichLat[i]:= TList.Create;

     // I have to go through all the cells within the desired latitude
     // Here I indicate the first and the last cell (first cell + total number of cells)
     For cell:=Trunc(CellLatMat[i,1]) to Trunc(CellLatMat[i,1]+CellLatMat[i,2])-1 do
     begin

       // For each cell I have to take all the species ID (one by one)
       For j:=0 to TRichCell[cell].Count-1 do
       begin

         // At each step I will take the species name and check if it is already in the species list of the latitude
         pSP:= TRichCell[cell][j];

         If TRichLat[i].IndexOf(pSP) < 0 then
         begin

           // If the species ID is not in the species list of the latitude, add it
           TRichLat[i].Add(pSP);
         end;

       end;

     end;

   end;

 end;

procedure equalSamplingApply(n: integer);
 var
 //sp: TIntVector;
 i, j, events, cell: Integer;

 begin

   // Vector with the number of species to be sampled at each sampling event (sp) was created as a global variable in this simulation

   // For each latitude I will randomly pick some cells to sample the species richness
   For i:=0 to Length(CellLatMat)-1 do
   begin

     // Start a TList for each latitude to save the ID of sampled species
     If Scont=0 then
     begin
       TObsLat[i]:= TList.Create;
     end;

     // The number of cells to pick is based on the number of replicates (parameter n)
     For events:=0 to n-1 do
     begin

       // The column 1 has the position of the first cell and the column 2 has the number of cells of each latitude
       // So I will draw only cells of the desired latitude
       // Note that it is allowed to sample the same cell many times, as occur in real world
       cell:= Trunc(CellLatMat[i,1])  + Random(Trunc(CellLatMat[i,2]) );

       // After draw a cell I will sample two species (two species ID) inside it (all cells have at least two species)
       // This number was based on the average number of species recorded by real sampling events
       sp[0]:= Random(TRichCell[cell].Count);
       sp[1]:= Random(TRichCell[cell].Count);
       // Sample again the second species in case of duplicate, but avoid a possible infinite loop (in case of a cell with only one species)
       If TRichCell[cell].Count>1 then
       begin
         While sp[1] = sp[0] do
         begin
           sp[1]:= Random(TRichCell[cell].Count);
         end;
       end;

       // For each sampled species...
       For j:=0 to Length(sp)-1 do
       begin

         // Save the species ID in the sampling cell TList (Database)
         pSP:= TRichCell[cell][sp[j]];
         TSampCell[cell].Add(pSP);

         // Save the species ID in the observed richness cell TList
         If TObsCell[cell].IndexOf(pSP) < 0 then
         begin
           TObsCell[cell].Add(pSP);
         end;

         // Save the species ID in the observed richness latitude TList
         If TObsLat[i].IndexOf(pSP) < 0 then
         begin
           TObsLat[i].Add(pSP);
         end;

       end;

     end;

   end;

 end;

procedure CellRichEstimates;
 var
 Sobs, Q1, Q2, m: Integer;
 i, j, k, contJ: Integer;

 begin

   // For each cell...
   For i:=0 to Length(neighMat)-1 do
   begin

     // Count the singletons and doubletons
     Q1:= 0;
     Q2:= 0;

     // First, pick one sampled species from the observed richness TList
     For j:=0 to TObsCell[i].Count-1 do
     begin

       // Then, count how many times this species was sampled at this cell (i.e., its incidence)
       contJ:= 0;
       For k:=0 to TSampCell[i].Count-1 do
       begin

         If TObsCell[i][j] = TSampCell[i][k] then
         begin
           contJ:= contJ + 1;
         end;

       end;

       // If it is a singleton or a doubleton, update the information
       If contJ = 1 then
       begin
         Q1:= Q1 + 1;
       end;
       If contJ = 2 then
       begin
         Q2:= Q2 + 1;
       end;

     end;

     // Calculate the estimated species richness
     Sobs:= TObsCell[i].Count;

     If Sobs>0 then
     begin

       If Q2=0 then
       begin

         // If there is no doubleton in this cell, calculate the Chao 2 'bias-corrected' version
         // Each sampling event colect two species. So, the number of samples is half the number of specimens
         m:= Trunc(TSampCell[i].Count/2);
         TEffort[Scont].MatEstCell0[i,nSimu]:= Sobs + ( ((m-1)/m) * ( (Q1*(Q1-1)) / (2*(Q2+1)) ) );

       end
       Else
       begin

         TEffort[Scont].MatEstCell0[i,nSimu]:= Sobs + (Power(Q1,2) / (2*Q2));

       end;

     end
     Else
     begin

       // If there is no observed species indicate with "-1", to differentiate of a richness estimate equal to zero
       TEffort[Scont].MatEstCell0[i,nSimu]:= -1;

     end;

   end;

 end;

procedure CoverageLat;
 var
 Sobs, Q1, Q2, m, n: Integer;
 i, j, k, cell, contJ: Integer;
 Q0_hat, A: Double;

 begin

   // This procedure calculate the inventory completeness based on sample coverage
   // But we will take to opportunity to estimate the especies richness by Chao 2 for the latitude

   // For each latitude...
   For i:=0 to Length(CellLatMat)-1 do
   begin

     // Count the singletons and doubletons
     Q1:= 0;
     Q2:= 0;

     // First, pick one sampled species from the observed richness TList
     For j:=0 to TObsLat[i].Count-1 do
     begin

       // Then, count how many times this species was sampled at this latitude (i.e., its incidence)
       contJ:= 0;
       For cell:=Trunc(CellLatMat[i,1]) to (Trunc(CellLatMat[i,1]) + Trunc(CellLatMat[i,2]))-1 do
       begin

         For k:=0 to TSampCell[cell].Count-1 do
         begin

           If TObsLat[i][j] = TSampCell[cell][k] then
           begin
             contJ:= contJ + 1;
           end;

         end;

       end;

       // If it is a singleton or a doubleton, update the information
       If contJ = 1 then
       begin
         Q1:= Q1 + 1;
       end;
       If contJ = 2 then
       begin
         Q2:= Q2 + 1;
       end;

     end;

     // Calculate the estimated species richness
     Sobs:= TObsLat[i].Count;

     If Sobs>0 then
     begin

       If Q2=0 then
       begin

         // If there is no doubleton in this latitude, calculate the Chao 2 'bias-corrected' version
         // The number of samples by latitude is already defined in the simulation
         m:= nSamp*(Scont+1);
         TEffort[Scont].MatEstLat0[i,nSimu]:= Sobs + ( ((m-1)/m) * ( (Q1*(Q1-1)) / (2*(Q2+1)) ) );

       end
       Else
       begin

         TEffort[Scont].MatEstLat0[i,nSimu]:= Sobs + (Power(Q1,2) / (2*Q2));

       end;

     end
     Else
     begin

       // If there is no observed species indicate with "-1", to differentiate of a richness estimate equal to zero
       TEffort[Scont].MatEstLat0[i,nSimu]:= -1;

     end;

     // Calculate completeness estimate (code inspired on iNEXT R package)

     // Define the number of records
     // Because each sampling event colect two species, the number of recors is the double of the number of sampling events
     n:= nSamp*(Scont+1) * 2;

     If Q2=0 then
     begin
       Q0_hat:= ((n-1)/n) * Q1 * ((Q1-1)/2);
     end
     Else
     begin
       Q0_hat:= ((n-1)/n) * Power(Q1,2)/2/Q2;
     end;

     If Q1>0 then
     begin
       A:= (n * Q0_hat)/((n * Q0_hat) + Q1);
     end
     Else
     begin
       A:= 1;
     end;

     TEffort[Scont].MatCoverageLat0[i,nSimu]:= 1 - ((Q1/n) * A);

   end;

 end;

procedure CountGapsLat;
 var
 i, j, MinPos, MaxPos: Integer;
 sum: Double;
 TUnic: TList;
 //vec: TIntVector;
 //MatGaps: TDblMatrix;

 begin

   // Vector to save presence/absence (vec) and the gaps matrix (MatGaps) were created as global variable in this simulation

   // Start a TList (details below)
   TUnic:= TList.Create;

   // The first step is to create TList with the ID of all sampled species (like a global checklist)
   // For each latitude...
   For i:=0 to Length(CellLatMat)-1 do
   begin

     // Check all species sampled
     For j:=0 to TObsLat[i].Count-1 do
     begin

       // If the species is not in the TList, add it
       pSP:= TObsLat[i][j];
       If TUnic.IndexOf(pSP) < 0 then
       begin
         TUnic.Add(pSP)
       end;

     end;

   end;

   // For each species...
   For i:=0 to TUnic.Count-1 do
   begin

     // Find where it occurs and define its range limits
     pSP:= TUnic[i];
     MinPos:= 26;
     MaxPos:= 0;

     // Run all the latitudes
     For j:=0 to Length(CellLatMat)-1 do
     begin

       If TObsLat[j].IndexOf(pSP) < 0 then
       begin

         // If the species was not sampled at this latitude, add 0
         vec[j]:= 0;
       end
       Else
       begin

         // Else, add 1
         vec[j]:= 1;

         // Update the species latitudinal range limits
         If j < MinPos then
         begin
           MinPos:= j;
         end;
         If j > MaxPos then
         begin
           MaxPos:= j;
         end;

       end;

     end;

     // Now, run the latitudinal range of the species to find the spatial gaps
     For j:=MinPos to MaxPos do
     begin
       If vec[j]=0 then
       begin
         MatGaps[j,i]:= 1;
       end;
     end;

   end;

   // After find the latitudinal gaps of all species, calculate the total number of gaps by latitude
   // For each latitude...
   For i:= 0 to Length(MatGaps)-1 do
   begin

     // Sum the number of species with spatial gaps (and clean the gaps matrix for the next sampling round)
     sum:= 0;
     For j:= 0 to Length(MatGaps[0])-1 do
     begin
       sum:= sum + MatGaps[i,j];
       MatGaps[i,j]:=0;
     end;

     TEffort[Scont].MatGapsLat0[i,nSimu]:= sum;
   end;

   // Finish the TList
   TUnic.Free

 end;

procedure rowSum(input: TDblMatrix; output: TDblMatrix; col: Integer);
 var
 sum: Double;
 i, j: Integer;

 begin

   // For each row of the input matrix (a latitude or cell)...
   For i:= 0 to Length(input)-1 do
   begin

     // Sum up the results of each simulation
     sum:= 0;
     For j:= 0 to Length(input[0])-1 do
     begin
       sum:= sum + input[i,j];
     end;

     // Calculate the mean and save in the output matrix
     // Each col in the output matrix is the mean result of a sampling scenario
     output[i,col]:= sum/Length(input[0]);
   end;

 end;


begin
  // Define input and output directory
  WriteLn('Set input directory path: ');
  ReadLn(inputDir);
  WriteLn('Set output directory path: ');
  ReadLn(outputDir);
  // Define the number of species and simulations
  WriteLn('Set number of species: ');
  ReadLn(Total_ssp);
  WriteLn('Set number of simulations: ');
  ReadLn(TotalSimu);
  // Define the range structure and maximum sampling effort
  WriteLn('Contiguous geographic range? y or n: ');
  ReadLn(Model);
  WriteLn('Set rate of increasing sampling effort: ');
  ReadLn(nSamp);
  WriteLn('Set maximum number of sampling events: ');
  ReadLn(maxSamp);

  // Input data
  SetCurrentDir(inputDir);
  ImportASCIIFile(neighMat, 'NeigCellMat.txt');
  ImportASCIIFile(CellLatMat, 'CellLatMat.txt');

  // Create ID of the cells (geographic space) and species (they are a reference to the pointers)
  // I use Cell_ID to indicate the cells within the range of a species, and SP_ID to indicate which species are within a cell
  SetLength(Cell_ID, Length(neighMat));
  For idx:=0 to Length(Cell_ID)-1 do
  begin
    Cell_ID[idx]:= idx;
  end;

  SetLength(SP_ID, Total_ssp);
  For idx:=0 to Total_ssp-1 do
  begin
    SP_ID[idx]:= idx;
  end;

  // Create array of TLists used along the simulation
  SetLength(TRichLat, Length(CellLatMat));
  SetLength(TObsLat, Length(CellLatMat));
  SetLength(TRichCell, Length(neighMat));
  SetLength(TObsCell, Length(neighMat));
  SetLength(TSampCell, Length(neighMat));

  // Create the generic variables of the equalSamplingApply and CountGapsLat procedures
  SetLength(sp, 2);
  SetLength(vec, Length(CellLatMat));
  SetLength(MatGaps, Length(CellLatMat), Total_ssp);

  // Create array of TSampling to save the information of each sampling scenario (the size depends on the number of scenarios)
  // This means that we have a matrix to save the result of each sampling scenario
  SetLength(TEffort, Trunc(maxSamp/nSamp));
  For idx:= 0 to Length(TEffort)-1 do
  begin
    SetLength(TEffort[idx].MatObsLat0, Length(CellLatMat), TotalSimu);
    SetLength(TEffort[idx].MatEstLat0, Length(CellLatMat), TotalSimu);
    SetLength(TEffort[idx].MatCoverageLat0, Length(CellLatMat), TotalSimu);
    SetLength(TEffort[idx].MatGapsLat0, Length(CellLatMat), TotalSimu);

    SetLength(TEffort[idx].MatObsCell0, Length(neighMat), TotalSimu);
    SetLength(TEffort[idx].MatEstCell0, Length(neighMat), TotalSimu);
  end;

  //Create matrices for output (In this version we'll save the mean of each sampling scenario)
  SetLength(MatRichLat, Length(CellLatMat), TotalSimu);
  SetLength(MatRichCell, Length(neighMat), TotalSimu);

  SetLength(MatObsLat, Length(CellLatMat), Length(TEffort));
  SetLength(MatEstLat, Length(CellLatMat), Length(TEffort));
  SetLength(MatCoverageLat, Length(CellLatMat), Length(TEffort));
  SetLength(MatGapsLat, Length(CellLatMat), Length(TEffort));
  SetLength(MatObsCell, Length(neighMat), Length(TEffort));
  SetLength(MatEstCell, Length(neighMat), Length(TEffort));


  // Start the simulation
  Randomize;

  For nSimu:= 0 to TotalSimu-1 do
  begin

    // Start a TList for each cell to save the ID of species occuping the cell (important to calculate spatial diversity)
    // Take the opportunity to start a similar TList to save the ID of species sampled
    For idx:=0 to Length(TRichCell)-1 do
    begin
      TRichCell[idx]:= TList.Create;
      TObsCell[idx]:= TList.Create;
      TSampCell[idx]:= TList.Create;
    end;

    // Run the range filling model for each species
    // You can choose a SpreadingDye (contiguous) or ScatteredRange model (not contiguous)
    // Both procedures run the model for one species and save the geographical distribution information
    If Model = 'y' then
    begin
      For Species:= 0 to Total_ssp-1 do
      begin
        SpreadingDye;
      end;
    end;
    If Model = 'n' then
    begin
      For Species:= 0 to Total_ssp-1 do
      begin
        ScatteredRange;
      end;
    end;

    // Calculate the known species richness of each latitude
    CountRichLat;

    // Save information (by latitude and cell) on output matrices (for known richness)
    For idx:=0 to Length(MatRichLat)-1 do
    begin
      MatRichLat[idx,nSimu]:= TRichLat[idx].Count;
    end;

    For idx:=0 to Length(MatRichCell)-1 do
    begin
      MatRichCell[idx,nSimu]:= TRichCell[idx].Count;
    end;


    // Simulate the sampling effort of each latitude
    For Scont:= 0 to Trunc(maxSamp/nSamp)-1 do
    begin

      // Apply the sampling effort simulation (Defined by the user)
      equalSamplingApply(nSamp);

      // Estimate the species richness using Chao 2, and calculate the inventory completeness (sample coverage)
      CellRichEstimates;
      CoverageLat;

      // Count the spatial gaps in species latitudinal distribution
      CountGapsLat;

      // Save information (by latitude and cell) on output matrices (estimates were already saved within the procedures)
      For idx:=0 to Length(MatObsLat)-1 do
      begin
        TEffort[Scont].MatObsLat0[idx,nSimu]:= TObsLat[idx].Count;
      end;

      For idx:=0 to Length(MatObsCell)-1 do
      begin
        TEffort[Scont].MatObsCell0[idx,nSimu]:= TObsCell[idx].Count;
      end;

    end;

    // Finish species TList by latitude and cell
    For idx:=0 to Length(TRichLat)-1 do
    begin
      TRichLat[idx].Free;
      TObsLat[idx].Free;
    end;

    For idx:=0 to Length(TRichCell)-1 do
    begin
      TRichCell[idx].Free;
      TObsCell[idx].Free;
      TSampCell[idx].Free;
    end;

  end;


  // Now calculate the mean result of each sampling scenario (for matrices of latitude and cell)
  For idx:= 0 to Length(TEffort)-1 do
  begin
    rowSum(TEffort[idx].MatObsLat0, MatObsLat, idx);
    rowSum(TEffort[idx].MatEstLat0, MatEstLat, idx);
    rowSum(TEffort[idx].MatCoverageLat0, MatCoverageLat, idx);
    rowSum(TEffort[idx].MatGapsLat0, MatGapsLat, idx);
  end;

  For idx:= 0 to Length(TEffort)-1 do
  begin
    rowSum(TEffort[idx].MatObsCell0, MatObsCell, idx);
    rowSum(TEffort[idx].MatEstCell0, MatEstCell, idx);
  end;

  // Export results
  SetCurrentDir(outputDir);
  QuickSaveData('MatRichLatES.txt',MatRichLat,Nil);
  QuickSaveData('MatObsLatES.txt',MatObsLat,Nil);
  QuickSaveData('MatEstLatES.txt',MatEstLat,Nil);
  QuickSaveData('MatCoverageLatES.txt',MatCoverageLat,Nil);
  QuickSaveData('MatGapsLatES.txt',MatGapsLat,Nil);

  QuickSaveData('MatRichCellES.txt',MatRichCell,Nil);
  QuickSaveData('MatObsCellES.txt',MatObsCell,Nil);
  QuickSaveData('MatEstCellES.txt',MatEstCell,Nil);

end.

