program SampSim;

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
  LibGammaFunc,
  LibProbDistrib;

Type
  Tsampling = Record
    MatObsLat0: TDblMatrix;
    MatCoverageLat0: TDblMatrix;
    MatChao2Lat0: TDblMatrix;
    MatInExtLat0: TDblMatrix;
    MatES50Lat0: TDblMatrix;
    MatGapsLat0: TDblMatrix;
    MatObsCell0: TDblMatrix;
    MatCoverageCell0: TDblMatrix;
    MatChao2Cell0: TDblMatrix;
    MatInExtCell0: TDblMatrix;
    MatES50Cell0: TDblMatrix;
   End;

  TCovPars = Record
    Q1: TIntVector;
    Q2: TIntVector;
    U: TIntVector;
   End;

var
  // Input matrices
  CellLatMat, neighMat: TDblMatrix;

  // User settings
  inputDir, outputDir: String;
  CohesiveRange, FounderLoc, SampType, RNDcells, recFlex: String;
  Total_ssp, TotalSimu, SampMin, SampMulti, fixRecs: Integer;
  aRSFD, bRSFD, aLoc, bLoc: Double;

  // Loop control (user settings / generic loops)
  Species, nSimu, Scont, idx: Integer;

  // TLists & pointers
  TRichLat, TObsLat, TRichCell, TObsCell: Array of TList;
  Cell_ID, SP_ID: TIntVector;
  pSP: ^Integer;

  //Temporary objects to save incidence frequency, sampling events and latitudinal gaps
  MatQLat, MatQCell, MatGaps: TDblMatrix;
  SEvLatVec, SEvCellVec, PresAbsVec: TIntVector;
  TCovParsLat, TCovParsCell: Array of TCovPars;

  // Output matrices
  TEffort: Array of Tsampling;
  MatRichLat, MatObsLat, MatCoverageLat, MatChao2Lat, MatInExtLat, MatES50Lat, MatGapsLat: TDblMatrix;
  MatRichCell, MatObsCell, MatCoverageCell, MatChao2Cell, MatInExtCell, MatES50Cell: TDblMatrix;
  MatObsCell2Lat, MatCoverageCell2Lat, MatChao2Cell2Lat, MatInExtCell2Lat, MatES50Cell2Lat: TDblMatrix;

const
  StatusMsg = #13'Processing %d of %d (%d%%) ...';

procedure SpreadingDye;
 var
 LimRange, Founder: Integer;
 i, x: Integer;
 TCell_Ocup, TCell_Neigh: TList;
 pOcup, pNeigh: ^Integer;

 begin

   // Find a random range size for this species (Note that the species must colonize at least 1 grid cell)
   LimRange:= Ceil(RndBeta(aRSFD, bRSFD)*Length(neighMat));
   If LimRange=0 then
   begin
     // Just for safety, because the function 'Ceil' should avoid a range = 0
     LimRange:=1;
   end;
   // Find a random grid cell to start the range filling
   Founder:= Floor(RndBeta(aLoc, bLoc)*Length(neighMat)); //Founder:= Random(Length(neighMat));

   // Save the ID of the species in the TList of this first grid cell (Remember to use 'SP_ID' as reference to the pointer because 'Species' change through the simulation)
   // This step is repeated below to each new grid cell occupied by the species
   pSP:=@SP_ID[Species];
   TRichCell[Founder].Add(pSP);

   // Create a TList of occupied and neighbors (i.e. available) grid cells
   TCell_Ocup:= TList.Create;
   TCell_Neigh:= TList.Create;
   // Save the starting point in the TList of occupied grid cells
   pOcup:=@Cell_ID[Founder];
   TCell_Ocup.Add(pOcup);

   // Start the range filling process
   While TCell_Ocup.Count <> LimRange do
   begin

     // First, find the neighbors grid cells
     // Here I will run all the rows of the neighborhood matrix to find the neighbors of the target grid cell (in the column)
     // The target grid cell is pointed by "pOcup", which is the Founder in the begining and randomly sampled after that
     // That is, I will always look for the neighbors of the last grid cell added to the list (i.e. occupied by the species)
     For i:= 0 to Length(neighMat)-1 do
     begin

       // If the grid cell is a neighbor, I will test if it can be occupied
       If (neighMat[i,pOcup^] = 1) then
       begin

         // It is important to use Cell_ID[] because it is my fixed reference
         // Note that I will use the memory address to test if the grid cell is already in the TLists
         pNeigh:= @Cell_ID[i];

         // Test if the grid cell (memory address) is already in the neighbors TList
         If TCell_Neigh.IndexOf(pNeigh) < 0 then
         begin

           // Test if the grid cell is already occupied
           If TCell_Ocup.IndexOf(pNeigh) < 0 then
           begin

             // If the neighbor cell is not occupied and it is not in the neighbors TList, add it to the neighbors TList
             // This way I will have a list only with potential grid cells for occupation
             TCell_Neigh.Add(pNeigh);
           end;

         end;

       end;

     end;

     // Now sample a neighbor cell to be occuppied
     x:= Random(TCell_Neigh.Count);

     // Again, I have to use the pointer because I'm working with memory address
     // Here I am transfering the address in the position 'x' to the pointer to add it in the list of occuppied grid cells
     pOcup:= TCell_Neigh[x];
     TCell_Ocup.Add(pOcup);

     // Delete the address of the now occupied grid cell from the neighbor TList
     TCell_Neigh.Delete(x);

     //Inform the presence of this species in the new occupied grid cell
     TRichCell[pOcup^].Add(pSP);
   end;

   // The species geographical range was created and the information is saved in the TRichCell TList
   // Finish the temporary TList
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

   // Find a randon range size for this species (Note that the species must colonize at least 1 grid cell)
   LimRange:= Ceil(RndBeta(aRSFD, bRSFD)*Length(neighMat)); //LimRange:= 1 + Random(Length(neighMat)-1);
   If LimRange=0 then
   begin
     // Just for safety, because the function 'Ceil' should avoid a range = 0
     LimRange:=1;
   end;
   // Find a randon grid cell to start the range filling
   Founder:= Floor(RndBeta(aLoc, bLoc)*Length(neighMat)); //Founder:= Random(Length(neighMat));

   // Save the ID of the species in the TList of this first grid cell (Remember to use 'SP_ID' as reference to the pointer because 'Species' change through the simulation)
   // This step is repeated below to each new grid cell occupied by the species
   pSP:=@SP_ID[Species];
   TRichCell[Founder].Add(pSP);

   // Create a TList of occupied grid cells (neighbors are not necessary here)
   TCell_Ocup:= TList.Create;
   // Save the starting point in the TList of occupied grid cells
   pOcup:=@Cell_ID[Founder];
   TCell_Ocup.Add(pOcup);

   // Start the range filling process
   While TCell_Ocup.Count <> LimRange do
   begin

     // Find a randon grid cell to continue the range filling
     x:= Random(Length(neighMat));
     pOcup:=@Cell_ID[x];

     // Test if the grid cell is already occupied
     If TCell_Ocup.IndexOf(pOcup) < 0 then
     begin

       // If the grid cell is not occupied, add it to the occupied grid cells TList
       TCell_Ocup.Add(pOcup);

       //Inform the presence of this species in the new occupied grid cell
       TRichCell[pOcup^].Add(pSP);
     end;

   end;

   // The species geographical range was created and the information is saved in the TRichCell TList
   // Finish the temporary TList
   TCell_Ocup.Free;
 end;

procedure CountRichness;
 var
 lat, cell, i: Integer;

 begin

   // For each latitude I will calculate the true species richness
   For lat:=0 to Length(CellLatMat)-1 do
   begin

     // Start a TList for each latitude to save the ID of all species present
     TRichLat[lat]:= TList.Create;

     // I have to go through all the grid cells within the focal latitude
     // Here I indicate the first and the last grid cell (first cell + total number of cells)
     For cell:=Trunc(CellLatMat[lat,1]) to Trunc(CellLatMat[lat,1]+CellLatMat[lat,2])-1 do
     begin

       // For each grid cell I have to take the ID of all species present (one by one)
       For i:=0 to TRichCell[cell].Count-1 do
       begin

         // At each step I will take the species ID and check if it is already in the species list of the latitude
         pSP:= TRichCell[cell][i];

         // If the species ID is not in the species list of the latitude, add it
         If TRichLat[lat].IndexOf(pSP) < 0 then
         begin
           TRichLat[lat].Add(pSP);
         end;

       end;

       // Take the opportunity to fill the MatRichCell output matrix
       MatRichCell[cell,nSimu]:= TRichCell[cell].Count;
     end;

     // Take the opportunity to fill the MatRichLat output matrix
     MatRichLat[lat,nSimu]:= TRichLat[lat].Count;
   end;

 end;

procedure SamplingApply(outMatLat, outMatCell: TDblMatrix);
 var
 TspSamp: TList;
 lat, cell, events, totEvents, nSP, sp, i: Integer;

 begin

   // Create a TList to save the ID of the species sampled at each sampling event
   TspSamp:= TList.Create;

   // For each latitude I will randomly pick some grid cells to sample the species
   For lat:=0 to Length(CellLatMat)-1 do
   begin

     // Start a TList for each latitude to save the ID of the sampled species
     If (Scont = 0) OR (SampType = 'y') then
     begin
       TObsLat[lat]:= TList.Create;
     end;

     // The number of grid cells to pick (with replacement) is based on the real number of sampling events at each latitudinal band
     // Or any sampling effort specified by the user when SampType = 'n'
     If SampType = 'y' then
     begin
       totEvents:= Trunc(CellLatMat[lat,3]);
     end
     Else
     begin

       // The first scenario executes the empirical sampling effort
       totEvents:= Trunc(CellLatMat[lat,3]);
       If (Scont > 0) then
       begin

         // After that, any effort below the minimal indicated by the user will be elevated
         // This minimal effort increases progressively at each sampling scenario
         If totEvents < (SampMin*Trunc(Power(2,Scont-1))) then
         begin
           totEvents:= (SampMin*Trunc(Power(2,Scont-1))) - SEvLatVec[lat];
         end
         Else
         begin
           totEvents:= 0;
         end;

       end;

     end;

     // At each sampling event...
     For events:=0 to totEvents-1 do
     begin

       // Draw a grid cell of the focal latitude
       // Column 1 has the position of the first grid cell and column 2 has the number of grid cells at each latitude
       // Note that it is allowed to sample the same grid cell multiple times, as occurs in the real world
       If RNDcells = 'y' then
       begin
         // This is the expected option
         cell:= Trunc(CellLatMat[lat,1]) + Random(Trunc(CellLatMat[lat,2]));
       end
       Else If RNDcells = 'n' then
       begin
          // This option was created only to explore the effects of the longitudinal distribution of sampling effort
          // This option only makes sense when there is a single row of cells per latitude
         cell:= Floor(RndBeta(CellLatMat[lat,7], CellLatMat[lat,7])*CellLatMat[lat,2]) + Trunc(CellLatMat[lat,1]);
       end;

       // Update the number of sampling events for the sampled grid cell (information necessary for Chao 2 richness estimation)
       SEvCellVec[cell]:= SEvCellVec[cell]+1;
       SEvLatVec[lat]:= SEvLatVec[lat]+1;

       // After draw a grid cell, I will draw a number of species (species ID) to sample inside it
       // This number is based on the distribution of species recorded by real sampling events at each latitude
       If recFlex = 'y' then
       begin
         nSP:= Ceil(RndBeta(CellLatMat[lat,4], CellLatMat[lat,5])*CellLatMat[lat,6]);
       end
       Else
       begin
         // Unless the user prefers to specify a fixed number of records (recFlex = 'n')
         nSP:= fixRecs;
       end;

       // If the number of species to sample is higher than the number of species present in the grid cell (unlikely), sample the species present.
       If TRichCell[cell].Count < nSP then
       begin
         nSP:= TRichCell[cell].Count;

         // If there is no species present in the grid cell (unlikely), go to the next sampling event
         If TRichCell[cell].Count = 0 then
         begin
           Continue;
         end;

       end;

       // Sampling the species
       While TspSamp.Count <> nSP do
       begin
         sp:= Random(TRichCell[cell].Count);
         pSP:= TRichCell[cell][sp];

         If TspSamp.IndexOf(pSP) < 0 then
         begin
           TspSamp.Add(pSP);
         end;

       end;

       // For each sampled species...
       For i:=0 to TspSamp.Count-1 do
       begin

         // Record the species ID in the specific "databases"
         pSP:= TspSamp[i];

         // Update the information in the incidence-frequency matrices
         MatQCell[cell,pSP^]:= MatQCell[cell,pSP^] + 1;
         MatQLat[lat,pSP^]:= MatQLat[lat,pSP^] + 1;

         // Save the species ID in the observed richness cell TList
         If TObsCell[cell].IndexOf(pSP) < 0 then
         begin
           TObsCell[cell].Add(pSP);
         end;

         // Save the species ID in the observed richness latitude TList
         If TObsLat[lat].IndexOf(pSP) < 0 then
         begin
           TObsLat[lat].Add(pSP);
         end;

       end;

       TspSamp.Clear;
     end;

     // Take the opportunity to fill the MatObsLat output matrix
     outMatLat[lat,nSimu]:= TObsLat[lat].Count;

     // Also fill the MatObsCell output matrix
     For cell:=Trunc(CellLatMat[lat,1]) to Trunc(CellLatMat[lat,1]+CellLatMat[lat,2])-1 do
     begin
       outMatCell[cell,nSimu]:= TObsCell[cell].Count;
     end;

   end;

   TspSamp.Free;
 end;

procedure CountGapsLat(outMat:TDblMatrix);
 var
 lat, sp, MinPos, MaxPos: Integer;
 sum: Double;
 TUnic: TList;

 begin

   // Start a TList (details below) to save presence/absence information
   TUnic:= TList.Create;

   // The first step is to fill the TList with the ID of all sampled species (like a global checklist)
   // For each latitude...
   For lat:=0 to Length(CellLatMat)-1 do
   begin

     // Check all species sampled
     For sp:=0 to TObsLat[lat].Count-1 do
     begin

       // If the species is not in the TList, add it
       pSP:= TObsLat[lat][sp];

       If TUnic.IndexOf(pSP) < 0 then
       begin
         TUnic.Add(pSP)
       end;

     end;

   end;

   // Now, for each sampled species...
   For sp:=0 to TUnic.Count-1 do
   begin

     // Find the latitudes where it occurs and define its range limits
     pSP:= TUnic[sp];
     MinPos:= Length(CellLatMat);
     MaxPos:= 0;

     // Run all the latitudes
     For lat:=0 to Length(CellLatMat)-1 do
     begin

       If TObsLat[lat].IndexOf(pSP) < 0 then
       begin
         // If the species was not sampled at this latitude, add 0 in the presence-absence vector
         PresAbsVec[lat]:= 0;
       end
       Else
       begin
         // Else, add 1
         PresAbsVec[lat]:= 1;

         // Update the species latitudinal range limits
         If lat < MinPos then
         begin
           MinPos:= lat;
         end;
         If lat > MaxPos then
         begin
           MaxPos:= lat;
         end;

       end;

     end;

     // Once we know its range limits, run the latitudinal range of the species to find spatial gaps
     For lat:=MinPos to MaxPos do
     begin
       If PresAbsVec[lat]=0 then
       begin
         MatGaps[lat,sp]:= 1;
       end;
     end;

   end;

   // After finding the latitudinal gaps of all species, calculate the total number of gaps at each latitude
   // For each latitude...
   For lat:= 0 to Length(MatGaps)-1 do
   begin

     // Sum the number of species with spatial gaps
     // Take the opportunity to reset the gaps matrix for the next sampling round
     sum:= 0;
     For sp:= 0 to Length(MatGaps[0])-1 do
     begin
       sum:= sum + MatGaps[lat,sp];
       MatGaps[lat,sp]:= 0;
     end;

     outMat[lat,nSimu]:= sum;
   end;

   // Finish the TList
   TUnic.Free
 end;

function CoverageFunc(inVec: TDblVector; Q1, Q2, U, T, tt: Integer): Double;
 var
 Q0_hat, A, sum, d: Double;
 sp: Integer;

 begin

   // Code inspired on iNEXT R package
   If Q2=0 then
   begin
     Q0_hat:= ((T-1)/T) * Q1 * ((Q1-1)/2);
   end
   Else
   begin
     Q0_hat:= ((T-1)/T) * Power(Q1,2)/2/Q2;
   end;

   If Q1>0 then
   begin
     A:= (T * Q0_hat)/((T * Q0_hat) + Q1);
   end
   Else
   begin
     A:= 1;
   end;

   If (T = tt) then
   begin
     // Conventional sample coverage estimate
     Result:= 1 - ((Q1/U) * A);
   end
   Else If (T > tt) then
   begin

     // Interpolation of the sample coverage estimate
     sum:= 0;
     For sp:=0 to Length(inVec)-1 do
     begin

       If (T - inVec[sp]) >= tt then
       begin
         sum:= sum + (inVec[sp]/U * exp(LnGamma(T-inVec[sp]+1, d)-LnGamma(T-inVec[sp]-tt+1, d)-LnGamma(T, d)+LnGamma(T-tt, d)));
       end;

     end;

     Result:= 1-sum;
   end
   Else If (T < tt) then
   begin
     // Extrapolation of the sample coverage estimate
     Result:= 1 - ((Q1/U) * Power(A,(tt-T+1)));
   end;

 end;

function Chao2Func(Q1, Q2, Sobs, T, tt: Integer): Double;
 var
 Q0_hat, A: Double;

 begin

   // If there is no doubleton in this reference sample, calculate the Chao 2 'bias-corrected' version
   If Q2=0 then
   begin
     Q0_hat:= ((T-1)/T) * Q1 * ((Q1-1)/2);
   end
   Else
   begin
     Q0_hat:= ((T-1)/T) * Power(Q1,2)/2/Q2;
   end;

   If Q1>0 then
   begin

     // This condition was included to account for extrapolation to a specific sampling size (tt)
     if tt <> 0 then
     begin
       // Originally the A equation is out of this condition, but when Q1=0 the value of A is undefined
       A:= (T * Q0_hat)/((T * Q0_hat) + Q1);
       Result:= Sobs + Q0_hat * (1 - Power(A, tt-T));
     end
     Else
     begin
       Result:= Sobs + Q0_hat;
     end;

   end
   Else
   begin
     Result:= Sobs;
   end;

 end;

function RarefyFunc(inVec: TDblVector; T, tt: integer): Double;
 var
 d, sum: Double;
 sp: Integer;

 begin

   // Calculate the estimated species richness
   If T >= tt then
   begin

     sum:= 0;
     For sp:=0 to Length(inVec)-1 do
     begin

       If (T-inVec[sp]) >= tt then
       begin
         // Note that the variable 'd' has no mathematical function here
         // It is only present because, for some reason, LnGamma demands a second double variable (with no function in the calculation)
         sum:= sum + (1 - exp(LnGamma(T-inVec[sp]+1, d)+LnGamma(T-tt+1, d)-LnGamma(T-inVec[sp]-tt+1, d)-LnGamma(T+1, d)));
       end
       Else
       begin
         // When T-inVec[sp] < tt, the resuting value for expression LnGamma(T-inVec[sp]-tt+1, d) will be zero or negative, resulting INF
         // In these cases the value can be turned to 0 (zero), and the result equals to 1.
         sum:= sum + 1;
       end;

     end;

     // Save the estimated richness in the output matrix
     Result:= sum;
   end
   Else
   begin
     Result:= -1; //NA
   end;

 end;

function FindSizeFunc(inVec: TDblVector; Q1, Q2, U, T, minT: Integer; refSC, stdSC:Double): Integer;
 var
 A, tempC, lowestC: Double;
 i, tt, contX, contInc, contDec: Integer;

 begin

   If refSC > stdSC then
   begin

     // In case of interpolation...
     // Find the sample size necessary to achieve the expected sample coverage

     // Shortcut to avoid unecessary time consuming runs
     // I assume here that the sampling effort necessary to achieve the expected coverage will be close to the sampling effort of the lowest coverage
     // But first I have to determine the direction
     lowestC:= 1;

     contDec:= 0;
     For i:= 0 to 2 do
     begin

       // First, lets check if the coverage becomes similar when the number of samples decreases
       tempC:= CoverageFunc(inVec, Q1, Q2, U, T, minT-i);
       tempC:= abs(tempC - stdSC);
       If tempC < lowestC then
       begin
         lowestC:= tempC;
         contDec:= contDec+1;
       end;

     end;

     contInc:= 0;
     For i:= 0 to 2 do
     begin

       // Then, lets check if the coverage becomes similar when the number of samples increases
       tempC:= CoverageFunc(inVec, Q1, Q2, U, T, minT+i);
       tempC:= abs(tempC - stdSC);
       If tempC < lowestC then
       begin
         lowestC:= tempC;
         contInc:= contInc+1;
       end;

     end;


     // Find the ideal sample size
     lowestC:= 1;
     contX:= 0;
     i:= minT;

     If contDec > contInc then
     begin

       // If decreasing sample size produces similar coverage more constantly...
       While contX = 0 do
       begin

         tempC:= CoverageFunc(inVec, Q1, Q2, U, T, i);
         tempC:= abs(tempC - stdSC);

         If (tempC < lowestC) AND (i <> 1) then
         begin
           lowestC:= tempC;
           tt:= i;
           i:= i-1;
         end
         Else
         begin
           // Once the tested coverage turns more dissimilar than similar to the reference, stop!
           contX:= contX + 1;
         end;

       end;

     end
     Else
     begin

       // Else, if increasing sample size produces similar coverage more constantly...
       While contX = 0 do
       begin

         tempC:= CoverageFunc(inVec, Q1, Q2, U, T, i);
         tempC:= abs(tempC - stdSC);

         If (tempC < lowestC) then
         begin
           lowestC:= tempC;
           tt:= i;
           i:= i+1;
         end
         Else
         begin
           // Once the tested coverage turns more dissimilar than similar to the reference, stop!
           contX:= contX + 1;
         end;

       end;

     end;

   end
   Else If refSC < stdSC then
   begin

     // In case of extrapolation...
     // Calculate the sample size necessary to achieve the expected sample coverage
     If (Q1>0) AND (Q2>0) then
     begin
       A:= ((T-1)*Q1) / ((T-1)*Q1 + 2*Q2);
     end
     Else If (Q1>1) AND (Q2=0) then
     begin
       A:= ((T-1)*(Q1-1)) / ((T-1)*(Q1-1) + 2);
     end
     Else
     begin
       A:= 1;
     end;

     tt:= Trunc( ((Ln(U/Q1) + Ln(1-stdSC)) / Ln(A)) - 1 );
     tt:= T+tt;
   end
   Else
   begin
     tt:= T;
   end;

   Result:= tt;
 end;

function InExtFunc(inVec: TDblVector; Q1, Q2, Sobs, T, tt: Integer): Double;
 begin

   If (T >= tt) then
   begin
     // In case of interpolation...
     Result:= RarefyFunc(inVec, T, tt);
   end
   Else
   begin
     // In case of extrapolation... apply Chao estimator for the specified sample size
     Result:= Chao2Func(Q1, Q2, Sobs, T, tt);
   end;

 end;

procedure EstimateDivCov(inMat: TDblMatrix; inTCovPars: Array of TCovPars; inObsVec: Array of TList; inSEvVec: TIntVector; outMatCov, outMatChao, outMatInExt, outMatES50: TDblMatrix);
 var
 i, sp, tt, minEffort: Integer;
 tempC, lowestC: Double;

 begin

   // For each sample reference, calculate the parameters necessary to execute the functions below
   For i:=0 to Length(inMat)-1 do
   begin

     inTCovPars[0].Q1[i]:= 0;
     inTCovPars[0].Q2[i]:= 0;
     inTCovPars[0].U[i]:= 0;

     // Count singletons and doubletons
     For sp:=0 to Length(inMat[0])-1 do
     begin

       If inMat[i,sp] = 1 then
       begin
         inTCovPars[0].Q1[i]:= inTCovPars[0].Q1[i] + 1;
       end;
       If inMat[i,sp] = 2 then
       begin
         inTCovPars[0].Q2[i]:= inTCovPars[0].Q2[i] + 1;
       end;

       // Number of records
       inTCovPars[0].U[i]:= inTCovPars[0].U[i] + Trunc(inMat[i,sp]);
     end;

   end;

   // Find the lowest sampling effort to interpolate richness for 'outMatES50' (unless it is lower than 50)
   minEffort:= 1000000000;
   For i:=0 to Length(inSEvVec)-1 do
   begin
     If (inSEvVec[i] > 0) AND (inSEvVec[i] <= minEffort) then
     begin
       minEffort:= inSEvVec[i];
       If minEffort < 50 then
       begin
         minEffort:= 50;
         break;
       end;
     end;
   end;

   // Calculate completeness and chao2 estimation
   For i:=0 to Length(inMat)-1 do
   begin

     // Condition suggested on Kusumoto et al (2020) to obtain reliable estimates
     If (inSEvVec[i]>5) AND (inTCovPars[0].Q1[i] <> inTCovPars[0].U[i]) then     //(inObsVec[i].Count>5) condition removed
     begin
       outMatCov[i,nSimu]:= CoverageFunc(inMat[i], inTCovPars[0].Q1[i], inTCovPars[0].Q2[i], inTCovPars[0].U[i], inSEvVec[i], inSEvVec[i]);
       outMatChao[i,nSimu]:= Chao2Func(inTCovPars[0].Q1[i], inTCovPars[0].Q2[i], inObsVec[i].Count, inSEvVec[i], 0);
       outMatES50[i,nSimu]:= RarefyFunc(inMat[i], inSEvVec[i], minEffort);
     end
     Else
     begin
       // If there was no sampling event or it was not possible to obtain reliable estimates, return "-1"
       outMatCov[i,nSimu]:= -1; //NA
       outMatChao[i,nSimu]:= -1; //NA
       outMatES50[i,nSimu]:= -1; //NA
     end;

   end;

   // Find the lowest sample coverage after doubling the sampling effort
   minEffort:= inSEvVec[0];
   lowestC:= 1;
   For i:=0 to Length(inMat)-1 do
   begin

     // Avoid locations where it was not possible to estimate sample coverage
     If outMatCov[i,nSimu] <> -1 then
     begin

       // If this sample coverage is the lowest found, save the information
       tempC:= CoverageFunc(inMat[i], inTCovPars[0].Q1[i], inTCovPars[0].Q2[i], inTCovPars[0].U[i], inSEvVec[i], 2*inSEvVec[i]);
       If tempC < lowestC then
       begin
         lowestC:= tempC;
         minEffort:= inSEvVec[i];
       end;

     end;

   end;

   // Once the optimal sample coverage reliable to species richness estimation was found, estimate species richness
   For i:=0 to Length(inMat)-1 do
   begin

     // Avoid locations where it was not possible to estimate sample coverage
     If outMatCov[i,nSimu] <> -1 then
     begin

       // Find the sample size necessary to achieve the expected sample coverage
       tt:= FindSizeFunc(inMat[i], inTCovPars[0].Q1[i], inTCovPars[0].Q2[i], inTCovPars[0].U[i], inSEvVec[i], minEffort, outMatCov[i,nSimu], lowestC);

       // Calculate the expected species richness for the standardized sample coverage
       outMatInExt[i,nSimu]:= InExtFunc(inMat[i], inTCovPars[0].Q1[i], inTCovPars[0].Q2[i], inObsVec[i].Count, inSEvVec[i], tt);
     end
     Else
     begin
       outMatInExt[i,nSimu]:= -1;
     end;

   end;

 end;

procedure CleanMat(inMat: TDblMatrix);
 var
 i, j: Integer;

 begin

   For i:= 0 to Length(inMat)-1 do
   begin

     For j:= 0 to Length(inMat[0])-1 do
     begin
       inMat[i,j]:= 0;
     end;

   end;

 end;

procedure rowMeans(inMat: TDblMatrix; outMat: TDblMatrix; col: Integer);
 var
 sum, mean, n: Double;
 i, j: Integer;

 begin

   // For each row of the input matrix (a latitude or grid cell)...
   For i:= 0 to Length(inMat)-1 do
   begin

     // Sum up the results of each simulation (the value in the columns)
     sum:= 0; n:= 0;
     For j:= 0 to Length(inMat[0])-1 do
     begin

       // Important to avoid cases when there was no sampling at the lat/cell and, therefore, species richness estimation is not available
       If inMat[i,j] > 0 then
       begin
         sum:= sum + inMat[i,j];
         n:= n + 1;
       end;

     end;

     // Calculate the mean value across the n simulations and save it in the output matrix
     // Each col in the output matrix is the mean, standard deviation and number of replicates of a specific sampling scenario
     If sum > 0 then
     begin

       // Number of replicates
       outMat[i,col+2]:= n;

       // Mean
       mean:= sum/n;
       outMat[i,col]:= mean;

       // Standard deviation
       If n > 1 then
       begin
         sum:= 0;
         For j:= 0 to Length(inMat[0])-1 do
         begin
           If inMat[i,j] > 0 then
           begin
             sum:= sum + Power((inMat[i,j]-mean),2);
           end;
         end;

         outMat[i,col+1]:= Sqrt(sum/(n-1));
       end
       Else
       begin
         outMat[i,col+1]:= -1;
       end;

     end
     Else
     begin
       outMat[i,col]:= -1;
       outMat[i,col+1]:= -1;
       outMat[i,col+2]:= -1;
     end;

   end;

 end;

procedure grid2lat(inMat: TDblMatrix; outMat: TDblMatrix; col: Integer);
 var
 sum, mean, n, sd, meanSum, num, den: Double;
 lat, simu, cell, meanCont, sdCont: Integer;

 begin
   // For each latitude...
   For lat:= 0 to Length(CellLatMat)-1 do
   begin

     meanSum:= 0; num:= 0; den:= 0; meanCont:= 0; sdCont:= 0;
     // At each simulation (the value in the columns)...
     For simu:= 0 to Length(inMat[0])-1 do
     begin

       // Sum up the results of the grid cells (the value in the rows)...
       sum:= 0; n:= 0;
       For cell:= Trunc(CellLatMat[lat,1]) to Trunc(CellLatMat[lat,1]+CellLatMat[lat,2])-1 do
       begin

         // Important to avoid cases when there was no sampling at the grid cell and, therefore, species richness estimation is not available
         If inMat[cell,simu] > 0 then
         begin
           sum:= sum + inMat[cell,simu];
           n:= n + 1;
         end;

       end;

       // Calculate the mean and SD for each latitude
       If sum > 0 then
       begin

         // Mean
         mean:= sum/n;
         meanSum:= meanSum + mean;
         meanCont:= meanCont + 1;

         // Standard deviation
         If n > 1 then
         begin
           sum:= 0;
           For cell:= Trunc(CellLatMat[lat,1]) to Trunc(CellLatMat[lat,1]+CellLatMat[lat,2])-1 do
           begin
             If inMat[cell,simu] > 0 then
             begin
               sum:= sum + Power((inMat[cell,simu]-mean),2);
             end;
           end;

           sd:= Sqrt(sum/(n-1));
           sdCont:= sdCont+1;

           num:= num + ( (n-1)*Power(sd, 2) );
           den:= den + n;
         end;

       end;

     end;

     // Calculate the averaged mean and averaged standard deviation
     If meanCont > 0 then
     begin
       outMat[lat,col]:= meanSum/meanCont;
     end
     Else
     begin
       outMat[lat,col]:= -1;
     end;

     If sdCont > 0 then
     begin
       den:= den - sdCont;
       outMat[lat,col+1]:= Sqrt(num/den);
       outMat[lat,col+2]:= sdCont;
     end
     Else
     begin
       outMat[lat,col+1]:= -1;
       outMat[lat,col+2]:= 0;
     end;

   end;

 end;

procedure ReportSettings;
 var
 arq: TextFile;

 begin
   AssignFile(arq, outputDir+'\Settings_Report.txt');
   Rewrite(arq);

   Write(arq, 'Number of species: ');
   WriteLn(arq, Total_ssp);

   Write(arq, 'Number of simulations: ');
   WriteLn(arq, TotalSimu);
   WriteLn(arq, '');
   WriteLn(arq, '');

   WriteLn(arq, '+--- Range Structure ---+');
   Write(arq, 'Contiguous geographic range?: ');
   WriteLn(arq, CohesiveRange);
   WriteLn(arq, '');

   Write(arq, 'Founder location random?: ');
   WriteLn(arq, FounderLoc);

   If FounderLoc = 'n' then
   begin
     Write(arq, #9 + 'Alpha and beta parameters to shape founder distribution: ');
     WriteLn(arq, FloatToStr(aLoc)+' '+FloatToStr(bLoc));
   end;

   WriteLn(arq, '');
   Write(arq, 'Alpha and beta parameters to shape the RSFD: ');
   WriteLn(arq, FloatToStr(aRSFD)+' '+FloatToStr(bRSFD));
   WriteLn(arq, '');
   WriteLn(arq, '');

   WriteLn(arq, '+--- Sampling Protocol ---+');
   Write(arq, 'Sampling effort latitudinally static?: ');
   WriteLn(arq, SampType);

   If SampType = 'n' then
   begin
     Write(arq, #9 + 'Minimal number of sampling events: ');
     WriteLn(arq, IntToStr(SampMin));
     Write(arq, #9 + 'Number of sampling scenarios: ');
     WriteLn(arq, IntToStr(SampMulti));
   end;

   WriteLn(arq, '');
   Write(arq, 'Records per sampling event flexible?: ');
   WriteLn(arq, recFlex);

   If recFlex = 'n' then
   begin
     Write(arq, #9 + 'Fixed number of records per sampling event: ');
     WriteLn(arq, IntToStr(fixRecs));
   end;

   WriteLn(arq, '');
   Write(arq, 'Sampling effort longitudinally random?: ');
   WriteLn(arq, RNDcells);

   CloseFile(arq);
 end;


begin
  // Define input and output directory
  WriteLn('Set input directory path: ');
  ReadLn(inputDir);
  WriteLn('');

  WriteLn('Set output directory path: ');
  ReadLn(outputDir);
  WriteLn('');

  // Define number of species and simulations
  WriteLn('Set number of species: ');
  ReadLn(Total_ssp);
  WriteLn('');

  WriteLn('Set number of simulations: ');
  ReadLn(TotalSimu);
  WriteLn('');
  WriteLn('');

  // Define range structure
  WriteLn('+--- Range Structure ---+');

  WriteLn('Contiguous geographic range? y or n: ');
  ReadLn(CohesiveRange);
  WriteLn('');

  WriteLn('Founder location random? y or n: ');
  ReadLn(FounderLoc);
  WriteLn('');

  If FounderLoc = 'n' then
  begin
    WriteLn(#9 + 'Set alpha and beta parameters to shape founder distribution: ');
    Write(#9); ReadLn(aLoc);
    Write(#9); ReadLn(bLoc);
    WriteLn('');
  end
  Else
  begin
    aLoc:= 1;
    bLoc:= 1;
  end;

  WriteLn('Set alpha and beta parameters to shape the Range Size Frequency Distribution: ');
  ReadLn(aRSFD);
  ReadLn(bRSFD);
  WriteLn('');
  WriteLn('');

  // Define sampling protocol
  WriteLn('+--- Sampling Protocol ---+');
  WriteLn('Sampling effort latitudinally static? y or n: ');
  ReadLn(SampType);
  WriteLn('');

  If SampType = 'n' then
  begin
    WriteLn(#9 + 'Set minimal number of sampling events: ');
    Write(#9); ReadLn(SampMin);
    WriteLn(#9 + 'Set number of sampling scenarios: ');
    Write(#9); ReadLn(SampMulti);
    WriteLn('');
  end;

  WriteLn('Records per sampling event flexible? y or n: ');
  ReadLn(recFlex);
  WriteLn('');

  If recFlex = 'n' then
  begin
    WriteLn(#9 + 'Set a fixed number of records per sampling event: ');
    Write(#9); ReadLn(fixRecs);
    WriteLn('');
  end;

  WriteLn('Sampling effort longitudinally random? y or n: ');
  ReadLn(RNDcells);
  WriteLn('');

  // Input data
  SetCurrentDir(inputDir);
  ImportASCIIFile(neighMat, 'NeigCellMat.txt');
  ImportASCIIFile(CellLatMat, 'CellLatMat.txt');

  // Create an ID for each grid cell (geographic space) and virtual species (the ID will be a reference to the pointers)
  // I use Cell_ID to indicate the grid cells within the range of a species, and SP_ID to indicate which species are within a grid cell
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

  // Create temporary matrices and vectors
  SetLength(MatQLat, Length(CellLatMat), Total_ssp);
  SetLength(MatQCell, Length(neighMat), Total_ssp);
  SetLength(MatGaps, Length(CellLatMat), Total_ssp);
  SetLength(SEvLatVec, Length(CellLatMat));
  SetLength(SEvCellVec, Length(neighMat));
  SetLength(PresAbsVec, Length(CellLatMat));

  SetLength(TCovParsLat, 1);
  For idx:= 0 to Length(TCovParsLat)-1 do
  begin
    SetLength(TCovParsLat[idx].Q1, Length(CellLatMat));
    SetLength(TCovParsLat[idx].Q2, Length(CellLatMat));
    SetLength(TCovParsLat[idx].U, Length(CellLatMat));
  end;

  SetLength(TCovParsCell, 1);
  For idx:= 0 to Length(TCovParsCell)-1 do
  begin
    SetLength(TCovParsCell[idx].Q1, Length(neighMat));
    SetLength(TCovParsCell[idx].Q2, Length(neighMat));
    SetLength(TCovParsCell[idx].U, Length(neighMat));
  end;

  // Create matrices for output
  SetLength(MatRichLat, Length(CellLatMat), TotalSimu);
  SetLength(MatRichCell, Length(neighMat), TotalSimu);

  If SampType = 'y' then
  begin
    // These matrices export the raw values of each simulation
    SetLength(MatObsLat, Length(CellLatMat), TotalSimu);
    SetLength(MatCoverageLat, Length(CellLatMat), TotalSimu);
    SetLength(MatChao2Lat, Length(CellLatMat), TotalSimu);
    SetLength(MatInExtLat, Length(CellLatMat), TotalSimu);
    SetLength(MatES50Lat, Length(CellLatMat), TotalSimu);
    SetLength(MatGapsLat, Length(CellLatMat), TotalSimu);

    SetLength(MatObsCell, Length(neighMat), TotalSimu);
    SetLength(MatCoverageCell, Length(neighMat), TotalSimu);
    SetLength(MatChao2Cell, Length(neighMat), TotalSimu);
    SetLength(MatInExtCell, Length(neighMat), TotalSimu);
    SetLength(MatES50Cell, Length(neighMat), TotalSimu);
  end
  Else If SampType = 'n' then
  begin

    // Create an array of TSampling to save the information of each sampling scenario
    // Each matrix in this array stores the results of one sampling scenario
    SetLength(TEffort, SampMulti);
    For idx:= 0 to Length(TEffort)-1 do
    begin
      SetLength(TEffort[idx].MatObsLat0, Length(CellLatMat), TotalSimu);
      SetLength(TEffort[idx].MatCoverageLat0, Length(CellLatMat), TotalSimu);
      SetLength(TEffort[idx].MatChao2Lat0, Length(CellLatMat), TotalSimu);
      SetLength(TEffort[idx].MatInExtLat0, Length(CellLatMat), TotalSimu);
      SetLength(TEffort[idx].MatES50Lat0, Length(CellLatMat), TotalSimu);
      SetLength(TEffort[idx].MatGapsLat0, Length(CellLatMat), TotalSimu);

      SetLength(TEffort[idx].MatObsCell0, Length(neighMat), TotalSimu);
      SetLength(TEffort[idx].MatCoverageCell0, Length(neighMat), TotalSimu);
      SetLength(TEffort[idx].MatChao2Cell0, Length(neighMat), TotalSimu);
      SetLength(TEffort[idx].MatInExtCell0, Length(neighMat), TotalSimu);
      SetLength(TEffort[idx].MatES50Cell0, Length(neighMat), TotalSimu);
    end;

    // These matrices export the average values across the simulations for each sampling scenario
    SetLength(MatObsLat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatCoverageLat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatChao2Lat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatInExtLat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatES50Lat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatGapsLat, Length(CellLatMat), Length(TEffort)*3);

    SetLength(MatObsCell, Length(neighMat), Length(TEffort)*3);
    SetLength(MatCoverageCell, Length(neighMat), Length(TEffort)*3);
    SetLength(MatChao2Cell, Length(neighMat), Length(TEffort)*3);
    SetLength(MatInExtCell, Length(neighMat), Length(TEffort)*3);
    SetLength(MatES50Cell, Length(neighMat), Length(TEffort)*3);

    SetLength(MatObsCell2Lat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatCoverageCell2Lat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatChao2Cell2Lat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatInExtCell2Lat, Length(CellLatMat), Length(TEffort)*3);
    SetLength(MatES50Cell2Lat, Length(CellLatMat), Length(TEffort)*3);
   end;


  // Start the simulation
  Randomize;

  For nSimu:= 0 to TotalSimu-1 do
  begin

    // Start a TList for each grid cell to save the ID of the species occuping the cells (necessary to run the range filling procedures)
    // Take the opportunity to start a similar TList to save the ID of the sampled species
    // Similar TLists corresponding to latitudes are started inside specific procedures (CountRichness & SamplingApply)
    For idx:=0 to Length(TRichCell)-1 do
    begin
      TRichCell[idx]:= TList.Create;
      TObsCell[idx]:= TList.Create;
    end;

    // Run the range filling model for each species
    // User can choose a SpreadingDye (contiguous) or ScatteredRange model (not contiguous)
    // Both procedures run the model for one species and save its geographical distribution information
    If CohesiveRange = 'y' then
    begin
      For Species:= 0 to Total_ssp-1 do
      begin
        SpreadingDye;
      end;
    end
    Else If CohesiveRange = 'n' then
    begin
      For Species:= 0 to Total_ssp-1 do
      begin
        ScatteredRange;
      end;
    end;

    // Calculate the true species richness for each grid cell and latitude (respective output matrices are called within the procedure)
    CountRichness;

    // Run the sampling simulation (and save the results in the respective output matrices)
    If SampType = 'y' then
    begin

      // Simulate the sampling effort for each latitude
      SamplingApply(MatObsLat, MatObsCell);

      // Count the spatial gaps in species latitudinal distribution
      CountGapsLat(MatGapsLat);

      // Calculate inventory completeness (sample coverage) and species richness estimation
      EstimateDivCov(MatQLat, TCovParsLat, TObsLat, SEvLatVec, MatCoverageLat, MatChao2Lat, MatInExtLat, MatES50Lat);
      EstimateDivCov(MatQCell, TCovParsCell, TObsCell, SEvCellVec, MatCoverageCell, MatChao2Cell, MatInExtCell, MatES50Cell);

    end
    Else If SampType = 'n' then
    begin

      // For each sampling scenario
      For Scont:= 0 to Length(TEffort)-1 do
      begin

        // Simulate the sampling effort for each latitude
        // Here, the incidence-frequency matrices are updated as sampling effort increases
        SamplingApply(TEffort[Scont].MatObsLat0, TEffort[Scont].MatObsCell0);

        // Count the spatial gaps in species latitudinal distribution
        CountGapsLat(TEffort[Scont].MatGapsLat0);

        // Calculate inventory completeness (sample coverage) and species richness estimation
        EstimateDivCov(MatQLat, TCovParsLat, TObsLat, SEvLatVec, TEffort[Scont].MatCoverageLat0, TEffort[Scont].MatChao2Lat0, TEffort[Scont].MatInExtLat0, TEffort[Scont].MatES50Lat0);
        EstimateDivCov(MatQCell, TCovParsCell, TObsCell, SEvCellVec, TEffort[Scont].MatCoverageCell0, TEffort[Scont].MatChao2Cell0, TEffort[Scont].MatInExtCell0, TEffort[Scont].MatES50Cell0);

      end;

    end;

    // Finish species TList by latitude and cell
    // Take the oportunity to reset temporary vectors of sampling effort
    For idx:=0 to Length(TRichLat)-1 do
    begin
      TRichLat[idx].Free;
      TObsLat[idx].Free;

      SEvLatVec[idx]:= 0;
    end;

    For idx:=0 to Length(TRichCell)-1 do
    begin
      TRichCell[idx].Free;
      TObsCell[idx].Free;

      SEvCellVec[idx]:= 0;
    end;

    // Reset temporary matrices (MatGaps & PresAbsVec are reset within the CountGapsLat procedure)
    CleanMat(MatQLat);
    CleanMat(MatQCell);

    // Progress bar
    Write(Format(StatusMsg, [nSimu+1, TotalSimu, Trunc(((nSimu+1)/TotalSimu)*100)]));
  end;


  // Export results
  If SampType = 'n' then
  begin

    // Calculate the mean result of each sampling scenario (for matrices of latitude and grid cells)
    Scont:= 0;
    For idx:= 0 to Length(TEffort)-1 do
    begin
      rowMeans(TEffort[idx].MatObsLat0, MatObsLat, Scont);
      rowMeans(TEffort[idx].MatCoverageLat0, MatCoverageLat, Scont);
      rowMeans(TEffort[idx].MatChao2Lat0, MatChao2Lat, Scont);
      rowMeans(TEffort[idx].MatInExtLat0, MatInExtLat, Scont);
      rowMeans(TEffort[idx].MatES50Lat0, MatES50Lat, Scont);
      rowMeans(TEffort[idx].MatGapsLat0, MatGapsLat, Scont);

      rowMeans(TEffort[idx].MatObsCell0, MatObsCell, Scont);
      rowMeans(TEffort[idx].MatCoverageCell0, MatCoverageCell, Scont);
      rowMeans(TEffort[idx].MatChao2Cell0, MatChao2Cell, Scont);
      rowMeans(TEffort[idx].MatInExtCell0, MatInExtCell, Scont);
      rowMeans(TEffort[idx].MatES50Cell0, MatES50Cell, Scont);

      grid2lat(TEffort[idx].MatObsCell0, MatObsCell2Lat, Scont);
      grid2lat(TEffort[idx].MatCoverageCell0, MatCoverageCell2Lat, Scont);
      grid2lat(TEffort[idx].MatChao2Cell0, MatChao2Cell2Lat, Scont);
      grid2lat(TEffort[idx].MatInExtCell0, MatInExtCell2Lat, Scont);
      grid2lat(TEffort[idx].MatES50Cell0, MatES50Cell2Lat, Scont);

      Scont:= Scont+3;
    end;

  end;

  SetCurrentDir(outputDir);
  QuickSaveData('MatRichLat'+'.txt',MatRichLat,Nil);
  QuickSaveData('MatObsLat'+'.txt',MatObsLat,Nil);
  QuickSaveData('MatCoverageLat'+'.txt',MatCoverageLat,Nil);
  QuickSaveData('MatChao2Lat'+'.txt',MatChao2Lat,Nil);
  QuickSaveData('MatInExtLat'+'.txt',MatInExtLat,Nil);
  QuickSaveData('MatES50Lat'+'.txt',MatES50Lat,Nil);
  QuickSaveData('MatGapsLat'+'.txt',MatGapsLat,Nil);

  QuickSaveData('MatRichCell'+'.txt',MatRichCell,Nil);
  QuickSaveData('MatObsCell'+'.txt',MatObsCell,Nil);
  QuickSaveData('MatCoverageCell'+'.txt',MatCoverageCell,Nil);
  QuickSaveData('MatChao2Cell'+'.txt',MatChao2Cell,Nil);
  QuickSaveData('MatInExtCell'+'.txt',MatInExtCell,Nil);
  QuickSaveData('MatES50Cell'+'.txt',MatES50Cell,Nil);

  If SampType = 'n' then
  begin
    QuickSaveData('MatObsCell2Lat'+'.txt',MatObsCell2Lat,Nil);
    QuickSaveData('MatCoverageCell2Lat'+'.txt',MatCoverageCell2Lat,Nil);
    QuickSaveData('MatChao2Cell2Lat'+'.txt',MatChao2Cell2Lat,Nil);
    QuickSaveData('MatInExtCell2Lat'+'.txt',MatInExtCell2Lat,Nil);
    QuickSaveData('MatES50Cell2Lat'+'.txt',MatES50Cell2Lat,Nil);
  end;

  ReportSettings;
end.


