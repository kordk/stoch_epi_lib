package main

import(
  "math"
  "math/rand"
  "sort"
	"os"
	"encoding/json"
	"io/ioutil"
  "time"
  "gonum.org/v1/gonum/mat"
  "strings"
  "fmt"
  "strconv"
  "flag"
  "github.com/montanaflynn/stats"
  "gonum.org/v1/gonum/integrate/quad"
  "runtime"
  "github.com/aclements/go-gg/generic/slice"
)


type Site struct{
  Lambda float64
  Lambdahat float64
  Mu float64
  Nu float64
  Alpha float64
  Kaploc []int
  Kapval []float64
  Aid string
}

type Gene struct{
  Gamma float64
  Dil float64
  Philoc []int
  Phival []float64
  GeneID string
}

type Realization struct{
  B [][]float64
  G [][]float64
  SG [][]float64
  T []float64
  GeneIDs []string
}

type MatchingData struct{
  Transcript map[string]float64
  Alpha map[string]float64
  Label string
}

type EqDist struct{
  GeneID string //gene ID
  DistEst []float64 //ith entry is f(x_i)
  XDomain []float64 // ith entry is x_i
  RealMeans []float64 // mean of the cumulative distribution at time point i
  RealVars []float64 // variance of the cumulative distribution at time point i
  TimePts []float64 // time points of cumulative mean/var
  MeanScaled float64
  StdScaled float64
}

type StateCount struct{
	State []float64
	TimeIn float64
}

type StateLikes struct{
  State []float64
  Li float64
}

type ParHelp2 struct{
  GrdSites []Site
  GrdGenes []Gene
  Likeli float64
  Label string
}


// HELPER FUNCTIONS TO MAKE COPIES.
// ------------------------------------------------------------------------------------

// Copy site parameters to new set of sites - used in TotalLike, StochGradDescent
func CopySites(oldsites []Site) []Site{
	newsites := make([]Site, len(oldsites))
	for i,s := range oldsites{
		newsites[i].Lambda = s.Lambda
		newsites[i].Lambdahat = s.Lambdahat
		newsites[i].Mu = s.Mu
		newsites[i].Alpha = s.Alpha
		newsites[i].Aid = s.Aid
		newsites[i].Kaploc = make([]int,len(s.Kaploc))
		copy(newsites[i].Kaploc,s.Kaploc)
		newsites[i].Kapval = make([]float64,len(s.Kapval))
		copy(newsites[i].Kapval,s.Kapval)
	}
	return newsites
}

// Copy gene parameters to new set of gene - used in TotalLike, StochGradDescent
func CopyGenes(oldgenes []Gene) []Gene{
	newgenes := make([]Gene, len(oldgenes))
	for i,g := range oldgenes{
		newgenes[i].Gamma = g.Gamma
		newgenes[i].Dil = g.Dil
		newgenes[i].GeneID = g.GeneID
		newgenes[i].Philoc = make([]int,len(g.Philoc))
		copy(newgenes[i].Philoc,g.Philoc)
		newgenes[i].Phival = make([]float64,len(g.Phival))
		copy(newgenes[i].Phival,g.Phival)
	}
	return newgenes
}

// Copy realization - used in Jumper, Jumper_use, Jumper2_use
func CopyRealization(r1 Realization) Realization{

	r2 := Realization{}
	r2.B = make([][]float64,len(r1.B))
	for i,x := range r1.B{
		r2.B[i] = make([]float64, len(x))
		copy(r2.B[i],x)
	}
	r2.G = make([][]float64,len(r1.G))
	for i,x := range r1.G{
		r2.G[i] = make([]float64, len(x))
		copy(r2.G[i],x)
	}
	r2.SG = make([][]float64,len(r1.SG))
	for i,x := range r1.SG{
		r2.SG[i] = make([]float64, len(x))
		copy(r2.SG[i],x)
	}
	r2.T = make([]float64,len(r1.T))
	copy(r2.T,r1.T)
	r2.GeneIDs = make([]string,len(r1.GeneIDs))
	copy(r2.GeneIDs,r1.GeneIDs)

	return r2
}


// Main helper funtions.
// ------------------------------------------------------------------------------------
//load model parameters from json file made in python
// Used in main
func BuildModel(sites_fl,genes_fl string) ([]Site, []Gene){
  	sitere, _ := os.Open(sites_fl)
  	siteby, _ := ioutil.ReadAll(sitere)
  	var mysiteli []Site
  	json.Unmarshal(siteby,&mysiteli)

  	genere, _ := os.Open(genes_fl)
  	geneby, _ := ioutil.ReadAll(genere)
  	var mygeneli []Gene
  	json.Unmarshal(geneby,&mygeneli)

    return mysiteli,mygeneli
}

// check if a file is in the working directory
// used in main
func FileExists(filename string) bool {
    info, err := os.Stat(filename)
    if os.IsNotExist(err) {
        return false
    }
    return !info.IsDir()
}



// PARAMETER FITTING FUNCTIONS
// ------------------------------------------------------------------------------------
//given a set of matching transcript/alpha values, try to optimize over possible parameter sets
//the log-likelihood of the data. Must be given an initial guess (which includes structural parameters it wont perturb)
//The alpha-values of the matching data should be given with the initial guess of sites
//This will try to optimize over Lambda,Lambdahat,Mu,Nu,Dil
// used in main
func MasterFit(sites []Site, genes []Gene, match_data []MatchingData, d ,bwidth float64, randoms [][]float64,ssteps int,stopnow int,numcrs int) ([]Site, []Gene, []float64) {
  fmt.Println("[MasterFit] Running.")
  //Which genes are going to behave "deterministically"?
  not_jumping := make([]int,0,len(genes))
  jumping_ones := make([]int,0,len(genes))
  mn_holder := make([]float64, len(match_data))
  all_means := make([]float64, len(genes))
  for i,gn := range genes{
    switch len(gn.Philoc){
    case 0:
      for j,dp := range match_data{
        mn_holder[j] = dp.Transcript[gn.GeneID]
      }
      meanof,_ := stats.Mean(mn_holder)
      all_means[i] = meanof
      not_jumping = append(not_jumping,i)
    default:
      for j,dp := range match_data{
        mn_holder[j] = dp.Transcript[gn.GeneID]
      }
      meanof,_ := stats.Mean(mn_holder)
      all_means[i] = meanof
      jumping_ones = append(jumping_ones,i)
    }
  }

  overall_mean,_ := stats.Mean(all_means)
  // overall_mean = overall_mean

  for _,i := range not_jumping{
    genes[i].Dil = genes[i].Dil*overall_mean
    genes[i].Gamma = all_means[i]*genes[i].Dil
  }

  gammabounds := make([]float64,len(jumping_ones))

  for l,i := range jumping_ones{
    genes[i].Dil = genes[i].Dil*overall_mean
    sm := 0.0
    for _,v := range genes[i].Phival{
      if v < 0{
        sm = sm - v
      }
    }
    gammabounds[l] = sm
  }

  // fmt.Println(gammabounds)

  fmt.Printf("[MasterFit] Fitting on %d/%d stochastic genes\n",len(jumping_ones),len(genes))
  //Now we have initialized on the surface with gn.Gamma/gn.Dil = mean(data[gn.GeneID]) for each of those.
  //We also don't want to change those parameters as we go, or solve for them or anything like that


  keepgoing := true // a stopping  condition here.
  num_nogood := 0
  fmt.Printf("[MasterFit] Computing initial likelihood.\n")

  first_like,numno,_ := TotalLike(sites,genes,match_data,randoms,bwidth,jumping_ones,overall_mean)
  fmt.Printf("[MasterFit] Initial likelihood: %0.5f\n No jump on %0.5f/%d draws.\n",first_like,numno,len(randoms[0]))
  likaswego := []float64{first_like}

  desc_tol := -6

  fmt.Println("[MasterFit] Generating sample parameter sets.")
  fmt.Println("[MasterFit] Stopping condition is %d steps in a row with descrease of less than 10^%d %",stopnow,desc_tol)

  mxstpszi := 0.2
  stpsize := 0.01
  for k:=0; (k<ssteps) && keepgoing; k++{
    fmt.Printf("[MasterFit] %d/%d Running. %s\n",k,ssteps, time.Now().Format("2006.01.02 15:04:05"))


    //we can do a data subset...that's the essence of ``stochastic" grad descent.
    subsize := 50
    var data_subset []MatchingData
    if len(match_data)<subsize{
      data_subset = match_data
    }else{
      data_subset = make([]MatchingData,subsize)
      r := rand.New(rand.NewSource(time.Now().Unix()))
      perm := r.Perm(len(match_data))
      for i := 0; i<subsize; i++{
        data_subset[i] = match_data[perm[i]]
      }
    }

    //do the grad descent
    new_sites,new_genes,gl2 := StochGradDescent(sites, genes, data_subset, d ,bwidth, randoms,stpsize,jumping_ones,overall_mean,numcrs)//gradient descent step needs to return parameters after 1 step
    if math.IsNaN(gl2){
      fmt.Println("[MasterFit] l2 norm of Gradient is now NaN, Stopping.")
      break
    }

    fmt.Printf("[MasterFit] found new parameters.\n")

    //compute the likelihood
    new_like,_,_ := TotalLike(new_sites,new_genes,match_data,randoms,bwidth,jumping_ones,overall_mean)//get the new likelihood

    fmt.Printf("[MasterFit] New likelihood %5f.\n", new_like)


    if math.IsNaN(new_like){
      fmt.Println("[MasterFit] New likelihood is now NaN, Stopping.")
      break
    }

    //as long as its not more I guess ok!
    if new_like < (1 - math.Pow(10,float64(desc_tol)))*likaswego[len(likaswego)-1]{
      num_nogood = 0
      sites = new_sites
      genes = new_genes

      if 2*stpsize < mxstpszi{
        stpsize = stpsize*1.1
      }
      likaswego = append(likaswego,new_like)
      fmt.Printf("[MasterFit] Likelihood now at %0.5f\n",new_like)

    }else{
      num_nogood += 1
      stpsize = stpsize*0.7
    }

    //quit if you haven't improved in 100 steps or whatever
    if num_nogood > stopnow{
      keepgoing = false
      fmt.Printf("[MasterFit] No improvement in %d steps, stopping", num_nogood)
    }
    if k+1 == ssteps{
      fmt.Printf("[MasterFit] Reached maximum allowed steps of %d",ssteps)
    }

  }

  for i,gn := range genes{
    genes[i].Dil = gn.Dil/overall_mean
  }

  // fmt.Println(likaswego)
  fmt.Printf("[MasterFit] Completed. %s\n",time.Now().Format("2006.01.02 15:04:05"))
  return sites,genes,likaswego
}

// Stochastic gradient descent of likelihood function of parameters.
// used in MasterFit
func StochGradDescent(sites []Site, genes []Gene, match_data []MatchingData, d ,bwidth float64, randoms [][]float64,tstep float64,includedlist []int,rescalefactor float64,numcrs int) ([]Site,[]Gene, float64){

  tmpsites := CopySites(sites)//make([]Site, 0)
  tmpgenes := CopyGenes(genes)//make([]Gene, 0)

  // numno := make([]float64,len(match_data))
  // totimes := make([]float64,len(match_data))

  bndwthmat := make([]float64,len(includedlist)*len(includedlist))
  for i := 0; i < len(includedlist); i++{
    for j := 0; j < len(includedlist); j++{
      if i==j{
        bndwthmat[i*len(includedlist) + j] = 1
      }
    }
  }

  topN := 100

  gradStDp := make(map[string][]Site, len(match_data))
  gradGnDp := make(map[string][]Gene, len(match_data))
  likesDp := make(map[string] float64, len(match_data))
  //
  switch numcrs {
  case 0:
    gradStDp,gradGnDp,likesDp = GradDescentDpSeries(sites, genes, tmpsites, tmpgenes, match_data, randoms,bwidth,bndwthmat, includedlist,rescalefactor, topN)
  case -1:
    nc := runtime.NumCPU()
    gradStDp,gradGnDp,likesDp = GradDescentDpParallel(sites, genes, tmpsites, tmpgenes, match_data, randoms,bwidth,bndwthmat, includedlist,rescalefactor, topN, nc)
  default:
    gradStDp,gradGnDp,likesDp = GradDescentDpParallel(sites, genes, tmpsites, tmpgenes, match_data, randoms,bwidth,bndwthmat, includedlist,rescalefactor, topN, numcrs)
  }


  gradVecSites := make([]Site,len(tmpsites))
  lngt := 0.0

  for j,_ := range tmpsites{
    //Lambda,Lambdahat,Mu,Nu
    minusGradLam := 0.0
    minusGradLamHat := 0.0
    minusGradMu := 0.0
    minusGradNu := 0.0

    for k,sitli := range gradStDp{
      minusGradLam += sitli[j].Lambda/likesDp[k]
      minusGradLamHat += sitli[j].Lambdahat/likesDp[k]
      minusGradMu += sitli[j].Mu/likesDp[k]
      minusGradNu += sitli[j].Nu/likesDp[k]
    }
    gradVecSites[j].Lambda = minusGradLam
    gradVecSites[j].Lambdahat = minusGradLamHat
    gradVecSites[j].Mu = minusGradMu
    gradVecSites[j].Nu = minusGradNu

    lngt += math.Pow(minusGradLam,2) + math.Pow(minusGradLamHat,2) + math.Pow(minusGradMu,2) + math.Pow(minusGradNu,2)
  }

  gradVecGns := make([]Gene,len(tmpgenes))

  for j,_ := range tmpgenes{
    //Dil
    minusGradDil := 0.0
    for k,gnli := range gradGnDp{
      minusGradDil += gnli[j].Dil/likesDp[k]
    }
    gradVecGns[j].Dil = minusGradDil
    lngt += math.Pow(minusGradDil,2)
  }

  Gradl2 := math.Sqrt(lngt)

  for j,si := range tmpsites{

    tmpsites[j].Lambda = math.Max(si.Lambda + tstep*gradVecSites[j].Lambda/Gradl2,0)
    tmpsites[j].Lambdahat = math.Max(si.Lambdahat + tstep*gradVecSites[j].Lambdahat/Gradl2,0)
    tmpsites[j].Mu = math.Max(si.Mu + tstep*gradVecSites[j].Mu/Gradl2,0)
    tmpsites[j].Nu = si.Nu + tstep*gradVecSites[j].Nu/Gradl2
  }


  for j,gn := range tmpgenes{
    tmpgenes[j].Dil = math.Max(gn.Dil + tstep*gradVecGns[j].Dil/Gradl2,0)
  }

  return tmpsites,tmpgenes,Gradl2

}

// Compute gradient for each datapoint in series.
// used in StochGradDescent
func GradDescentDpSeries(sites []Site, genes []Gene, tmpsites []Site, tmpgenes []Gene, match_data []MatchingData, randoms [][]float64,bwidth float64,bndwthmat []float64, includedlist []int,rescalefactor float64, topN int)(map[string][]Site,map[string][]Gene,map[string] float64){

  gradStDp := make(map[string][]Site, len(match_data))
  gradGnDp := make(map[string][]Gene, len(match_data))
  likesDp := make(map[string] float64, len(match_data))

  for _,dp := range match_data{
    gradStDp[dp.Label],gradGnDp[dp.Label],likesDp[dp.Label] = GradDescentDp(sites,genes,tmpsites,tmpgenes,dp,randoms,bwidth,bndwthmat,includedlist,rescalefactor, topN)
  }

  return gradStDp,gradGnDp,likesDp
}

// Compute gradient for each datapoint in parallel
// used in StochGradDescent
func GradDescentDpParallel(sites []Site, genes []Gene, tmpsites []Site, tmpgenes []Gene, match_data []MatchingData, randoms [][]float64,bwidth float64,bndwthmat []float64, includedlist []int,rescalefactor float64, topN int, numcores int)(map[string][]Site,map[string][]Gene,map[string] float64){


  runtime.GOMAXPROCS(numcores)


  gradStDp := make(map[string][]Site, len(match_data))
  gradGnDp := make(map[string][]Gene, len(match_data))
  likesDp := make(map[string] float64, len(match_data))


  i:= 0

  var ch chan ParHelp2

  for i < len(match_data){
    todo := int(math.Min(float64(len(match_data) - i),float64(numcores)))//how many to run at once.
    ch = make(chan ParHelp2,todo)

    for j:=0; j<todo; j++{
      go GradDescentParWrap(sites, genes, tmpsites, tmpgenes, match_data[i+j], randoms, bwidth, bndwthmat, includedlist, rescalefactor, topN, ch)
    }

    for jj := 0; jj < todo; jj++{
      var chv ParHelp2
      chv = <- ch
      gradStDp[chv.Label] = chv.GrdSites
      gradGnDp[chv.Label] = chv.GrdGenes
      likesDp[chv.Label] = chv.Likeli
    }
    i += todo
  }


  return gradStDp,gradGnDp,likesDp
}

// Wrapper for GradDescentDp so that it can be run in parallel
// used in GradDescentDpParallel
func GradDescentParWrap(sites []Site, genes []Gene, tmpsites []Site, tmpgenes []Gene, dp MatchingData, randoms [][]float64,bwidth float64,bndwthmat []float64, includedlist []int,rescalefactor float64, topN int, ch chan <- ParHelp2){


  gS,gG,L := GradDescentDp(sites, genes, tmpsites, tmpgenes, dp, randoms, bwidth, bndwthmat, includedlist, rescalefactor, topN)
  retval := ParHelp2{gS,gG,L,dp.Label}

  ch <- retval
}

// Compute gradient for a single data point
// used in GradDescentDpSeries, GradDescentParWrap
func GradDescentDp(sites []Site, genes []Gene, tmpsites []Site, tmpgenes []Gene, dp MatchingData, randoms [][]float64,bwidth float64,bndwthmat []float64, includedlist []int,rescalefactor float64, topN int)([]Site,[]Gene,float64){

  for j := 0; j< len(tmpsites); j++{
    tmpsites[j].Alpha = dp.Alpha[tmpsites[j].Aid]
  }
  //that 1000 parameter - maybe could be lower.
  stateCounts,StateLikelihoods := PrepStates(sites,genes,dp.Transcript,randoms,bwidth,bndwthmat,includedlist,rescalefactor, topN)

  likelihd := 0.0
  for _,B := range StateLikelihoods{
    likelihd += B.Li
  }

  transc_val := make([]float64,len(genes))
  for ii,gnn := range genes{
    transc_val[ii] = dp.Transcript[gnn.GeneID]
  }

  gradSt := make([]Site,len(tmpsites))
  gradGn := make([]Gene,len(tmpgenes))
  for j,si := range tmpsites{
    gradSt[j].Lambda = ComputeDerivSite(1, j,stateCounts, si,StateLikelihoods,transc_val)
    gradSt[j].Lambdahat = ComputeDerivSite(2, j,stateCounts, si,StateLikelihoods,transc_val)
    gradSt[j].Mu = ComputeDerivSite(3, j,stateCounts, si,StateLikelihoods,transc_val)
    gradSt[j].Nu = ComputeDerivSite(4, j,stateCounts, si,StateLikelihoods,transc_val)
  }
  for j,gn := range tmpgenes{
    gradGn[j].Dil = ComputeDerivGene(j,stateCounts,gn,dp.Transcript[gn.GeneID])
  }


  return gradSt,gradGn,likelihd

}

// Compute derivative for a gene parameter.
// used in GradDescentDp
func ComputeDerivGene(paramindex int, topStates []StateCount, gene Gene, transcriptval float64) float64{
  weighted_avg_deriv := 0.0
  for _,v := range topStates{

    //optimal dn = (1/gn)*(gamma_n + phiN dot B)
    phiDotB := 0.0
    for j,l := range gene.Philoc{
      phiDotB += gene.Phival[j]*v.State[l]
    }

    var optD float64
    if transcriptval > 0{
      optD = (1/transcriptval)*(gene.Gamma + phiDotB)
    }else{
      optD = 2*gene.Dil
    }

    diff := optD - gene.Dil

    weighted_avg_deriv += v.TimeIn*diff
  }
  return weighted_avg_deriv
}

// Compute derivative for a site parameter.
// used in GradDescentDp
func ComputeDerivSite(paramtype int, paramindex int,topStates []StateCount, site Site,all_partlikes []StateLikes, transcriptval []float64) float64 {
  //paramtypes are 1:Lambda, 2:Lambdahat, 3:Mu, 4:Nu
  // rhs := make([]float64,topN)

  topStatesSli := make([][]float64,len(topStates))
  for k,tS := range topStates{
    topStatesSli[k] = tS.State
  }

  sumDeriv := 0.0
  doneAlready := make([]int,0,len(topStates))

  gdotKap := 0.0
  for kvind,siind := range site.Kaploc{
    gdotKap += math.Max(transcriptval[siind],0.01)*site.Kapval[kvind]//we can actually use any value of g...as long as it isn't 0!
  }





  paramcomboFull := site.Lambda*gdotKap*site.Mu/(site.Mu + math.Pow(site.Alpha,site.Nu))

  switch paramtype{
  case 1://lambda
    paramcombo := gdotKap*site.Mu/(site.Mu + math.Pow(site.Alpha,site.Nu))
    for j,B := range topStatesSli{

      if NotIn(j,doneAlready){

        offb1 := GetOffByOne(B,topStatesSli, paramindex)//index of j Delta m

        if offb1 != -1{
          Bd := topStatesSli[offb1]

          // 2x2 system is
          // [A1 B1] [dl1dth] = [C1]
          // [A2 B2] [dldth]    [C2]

          C1 := paramcombo*B[paramindex]*(1-Bd[paramindex])*all_partlikes[offb1].Li - paramcombo*(1-B[paramindex])*all_partlikes[j].Li
          C2 := paramcombo*Bd[paramindex]*(1-B[paramindex])*all_partlikes[j].Li - paramcombo*(1-Bd[paramindex])*all_partlikes[offb1].Li

          A1 := -(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]
          A2 := -(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]

          B1 := B[paramindex]*(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]*(1-B[paramindex])
          B2 := Bd[paramindex]*(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]*(1-Bd[paramindex])



          sumDeriv += (1/(A1*B2 - A2*B1))*(B2*C1 - B1*C2 -A2*C1 + A1*C2)

          doneAlready = append(doneAlready,j)
          doneAlready = append(doneAlready,offb1)

        }else{

          C1 := - paramcombo*(1-B[paramindex])*all_partlikes[j].Li

          A1 := -(1-B[paramindex])*site.Lambda*paramcombo + site.Lambdahat*B[paramindex]

          sumDeriv += C1/A1

          doneAlready = append(doneAlready,j)



        }
      }
    }
  case 2://lambdahat

    for j,B := range topStatesSli{

      if NotIn(j,doneAlready){

        offb1 := GetOffByOne(B,topStatesSli, paramindex)//index of j Delta m

        if offb1 != -1{
          Bd := topStatesSli[offb1]

          // 2x2 system is
          // [A1 B1] [dl1dth] = [C1]
          // [A2 B2] [dldth]    [C2]

          C1 := Bd[paramindex]*(1-B[paramindex])*all_partlikes[offb1].Li - B[paramindex]*all_partlikes[j].Li
          C2 := B[paramindex]*(1-Bd[paramindex])*all_partlikes[j].Li - Bd[paramindex]*all_partlikes[offb1].Li

          A1 := -(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]
          A2 := -(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]

          B1 := B[paramindex]*(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]*(1-B[paramindex])
          B2 := Bd[paramindex]*(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]*(1-Bd[paramindex])



          sumDeriv += (1/(A1*B2 - A2*B1))*(B2*C1 - B1*C2 -A2*C1 + A1*C2)

          doneAlready = append(doneAlready,j)
          doneAlready = append(doneAlready,offb1)

        }else{

          C1 := -  B[paramindex]*all_partlikes[j].Li

          A1 := -(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]

          sumDeriv += C1/A1

          doneAlready = append(doneAlready,j)



        }
      }
    }

<<<<<<< HEAD
=======
    //as long as its not more I guess ok!
    if new_like < (1 - math.Pow(10,float64(desc_tol)))*likaswego[len(likaswego)-1]{
      num_nogood = 0
      sites = new_sites
      genes = new_genes
>>>>>>> d3a2cf5589242c4d24cb642dc30fd3d286bea85a


  case 3://mu
    paramcombo := gdotKap*(math.Pow(site.Alpha,site.Nu))/(site.Mu + math.Pow(site.Alpha,site.Nu))


    for j,B := range topStatesSli{

      if NotIn(j,doneAlready){

        offb1 := GetOffByOne(B,topStatesSli, paramindex)//index of j Delta m

        if offb1 != -1{
          Bd := topStatesSli[offb1]

          // 2x2 system is
          // [A1 B1] [dl1dth] = [C1]
          // [A2 B2] [dldth]    [C2]

          C1 := paramcombo*B[paramindex]*(1-Bd[paramindex])*all_partlikes[offb1].Li - paramcombo*(1-B[paramindex])*all_partlikes[j].Li
          C2 := paramcombo*Bd[paramindex]*(1-B[paramindex])*all_partlikes[j].Li - paramcombo*(1-Bd[paramindex])*all_partlikes[offb1].Li

          A1 := -(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]
          A2 := -(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]

          B1 := B[paramindex]*(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]*(1-B[paramindex])
          B2 := Bd[paramindex]*(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]*(1-Bd[paramindex])



          sumDeriv += (1/(A1*B2 - A2*B1))*(B2*C1 - B1*C2 -A2*C1 + A1*C2)

          doneAlready = append(doneAlready,j)
          doneAlready = append(doneAlready,offb1)

        }else{

          C1 := - paramcombo*(1-B[paramindex])*all_partlikes[j].Li

          A1 := -(1-B[paramindex])*site.Lambda*paramcombo + site.Lambdahat*B[paramindex]

          sumDeriv += C1/A1

          doneAlready = append(doneAlready,j)



        }
      }
    }



  case 4:

    sAlpha := math.Max(0.001,site.Alpha)

    paramcombo := -gdotKap*(site.Mu*math.Pow(sAlpha,site.Nu)*math.Log(sAlpha))/(math.Pow(site.Mu + math.Pow(sAlpha,site.Nu),2))



    for j,B := range topStatesSli{

      if NotIn(j,doneAlready){

        offb1 := GetOffByOne(B,topStatesSli, paramindex)//index of j Delta m

        if offb1 != -1{
          Bd := topStatesSli[offb1]

          // 2x2 system is
          // [A1 B1] [dl1dth] = [C1]
          // [A2 B2] [dldth]    [C2]

          C1 := paramcombo*B[paramindex]*(1-Bd[paramindex])*all_partlikes[offb1].Li - paramcombo*(1-B[paramindex])*all_partlikes[j].Li
          C2 := paramcombo*Bd[paramindex]*(1-B[paramindex])*all_partlikes[j].Li - paramcombo*(1-Bd[paramindex])*all_partlikes[offb1].Li

          A1 := -(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]
          A2 := -(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]

          B1 := B[paramindex]*(1-Bd[paramindex])*paramcomboFull + site.Lambdahat*Bd[paramindex]*(1-B[paramindex])
          B2 := Bd[paramindex]*(1-B[paramindex])*paramcomboFull + site.Lambdahat*B[paramindex]*(1-Bd[paramindex])



          sumDeriv += (1/(A1*B2 - A2*B1))*(B2*C1 - B1*C2 -A2*C1 + A1*C2)

          doneAlready = append(doneAlready,j)
          doneAlready = append(doneAlready,offb1)

        }else{

          C1 := - paramcombo*(1-B[paramindex])*all_partlikes[j].Li

          A1 := -(1-B[paramindex])*site.Lambda*paramcombo + site.Lambdahat*B[paramindex]

          sumDeriv += C1/A1

          doneAlready = append(doneAlready,j)



        }
      }
    }

  }
  return sumDeriv

}

// Run a sim and get data - top states, likelihoods by state
// used in GradDescentDp
func PrepStates(sites []Site, genes []Gene,transcript map[string]float64, randoms [][]float64,bndwthsiz float64, bndwthmat []float64, includedlist []int,rescalefactor float64, topN int) ([]StateCount,[]StateLikes) {
  //Need to get our realization into equilibrium
  //first need some initial conditions
  num_sites := len(sites)
  num_genes := len(genes)
  binit := make([]float64,num_sites)
  ginit := make([]float64,num_genes)
  transc_vals := make([]float64,len(includedlist))
  for i,gn := range genes{
    ginit[i] = transcript[gn.GeneID]/rescalefactor
  }
  for i,sj := range includedlist{
    transc_vals[i] = transcript[genes[sj].GeneID]/rescalefactor//ginit[sj]
  }
  //then get a sampler going.
  path,numnoj := Jumper_use(binit,ginit,0,sites,genes,randoms,"none")
  if float64(numnoj)/float64(len(randoms)) > 0.75{
    path,_ = Jumper2_use(binit,ginit,0,sites,genes,randoms,"none")
  }

  //Determine which site-states to use in solving the system
  states := make([]StateCount,0,len(path.T))
  for i:=0; i < len(path.T)-1;i++{
    dt := path.T[i+1]-path.T[i]
    state := path.B[i]
    appnd := true
    for i,v := range states{
			if SameSlice(state,v.State){
				states[i].TimeIn = states[i].TimeIn + dt
				appnd = false
			}
		}
		if appnd{
			states = append(states,StateCount{state,dt})
		}
  }


  topStates := TopNStates(states,topN,0.00001)


  topStatesSli := make([][]float64,len(topStates))
  for k,tS := range topStates{
    topStatesSli[k] = tS.State
  }

  //now get all the Li, dLidg for those
  all_partlikes := make([]StateLikes, len(topStates))
  for i,st := range topStates{
    all_partlikes[i] = GetLiState(path,st.State,genes,transc_vals,bndwthsiz,bndwthmat,includedlist)
  }
  return topStates,all_partlikes
}

// Get likelihood for single state.
// Used in PrepStates
func GetLiState(path Realization,Bstate []float64, genes []Gene,datapoint []float64,bndwthsiz float64, bndwthmat []float64, includedlist []int)(StateLikes){
  estim := 0.0

  qdpts := 3

  totalT := path.T[len(path.T)-1] - path.T[1]
  all_dils := make([]float64, len(includedlist))
  for i,sj := range includedlist{
    all_dils[i] = genes[sj].Dil
  }
  // fmt.Println(all_dils)
  all_giks := make([]float64,len(includedlist))
  all_siks := make([]float64,len(includedlist))
  for j:=1; j<len(path.T)-1; j++{
    if SameSlice(path.B[j],Bstate){
      Dtk := path.T[j+1]-path.T[j]
      // copy(all_giks,path.G[j])
      for i,sj := range includedlist{
        all_giks[i] = path.G[j][sj]
      }
      // copy(all_siks,path.SG[j])
      for i,sj := range includedlist{
        all_siks[i] = path.SG[j][sj]
      }
      eval := GaussKernal(all_dils,all_giks,all_siks,datapoint,bndwthmat,bndwthsiz)
      switch{
      case Dtk > 0:
        estim += quad.Fixed(eval, 0, Dtk, qdpts, nil, 0)
      case Dtk == 0:
        estim += 0
      case Dtk<0:
        panic("Time is moving backwards!!!")
      }
    }
  }
  return StateLikes{Bstate,(1.0/totalT)*estim}
}

// Checks if slices of floats are identical or not
// used in PrepStates, GetLiState
func SameSlice(x []float64,y []float64) bool{
	same := true
	if len(x) == len(y){
		for i:=0; i< len(x); i++{
			if x[i] != y[i]{
				same = false
			}
		}
	}else{
		same = false
	}
	return same
}

// Gets top N visted states by time
// used in PrepStates
func TopNStates(states []StateCount, N int, minTime float64) []StateCount{

	topN := make([]StateCount,0,len(states))
	just_Ts := make([]float64,0,len(states))

  ct := 0

  if N > 0{
    for len(topN) < N && ct < len(states){
      if states[ct].TimeIn > minTime{
        topN = append(topN,states[ct])
        just_Ts = append(just_Ts,states[ct].TimeIn)
      }
      ct += 1
    }

  	for i := N; i<len(states); i++{
  		aMin := slice.ArgMin(just_Ts)
  		if states[i].TimeIn > just_Ts[aMin]{
  			topN[aMin] = states[i]
  			just_Ts[aMin] = states[i].TimeIn
  		}
  	}
  }else{
    for ct < len(states){
      if states[ct].TimeIn > minTime{
        topN = append(topN,states[ct])
        just_Ts = append(just_Ts,states[ct].TimeIn)
      }
      ct += 1
    }
  }

	return topN
}

// Get site thats only different in index diffOn from Bst
// used in ComputeDerivSite
func GetOffByOne(Bst []float64, all_Bstates [][]float64, diffOn int) int {
  the_one := -1
  for j,bs := range all_Bstates{
    all_others_same := true
    this_one_diff := true
    for i,v := range bs{
      switch i{
      case diffOn:
        if v == Bst[i]{
          this_one_diff = false
        }
      default:
        if v != Bst[i]{
          all_others_same = false
        }
      }
    }
    if all_others_same && this_one_diff{
      the_one = j
      break
    }
  }
  return the_one
}

// Check if int a is in slice B
// used in ComputeDerivSite
func NotIn(a int, B []int) bool{
  NotIn := true
  for _,j := range B{
    if a == j{
      NotIn = false
      break
    }
  }
  return NotIn
}


// LIKELIHOOD COMPUTING FUNCTIONS
// ------------------------------------------------------------------------------------

// Get Likelihood of data according to a realization
// Used in MasterFit, main
func TotalLike(sites []Site, genes []Gene, data_points []MatchingData, randoms [][]float64,bndwthsiz float64,includedlist []int,rescalefactor float64) (float64,float64,float64){
  all_likes := make([]float64,len(data_points))
  tmpsites := CopySites(sites)//make([]Site, 0)

  //copier.Copy(&tmpsites,&sites)
  tmpgenes := CopyGenes(genes)//make([]Gene, 0)
  //copier.Copy(&tmpgenes,&genes)
  numno := make([]float64,len(data_points))
  totimes := make([]float64,len(data_points))

  bndwthmat := make([]float64,len(includedlist)*len(includedlist))
  for i := 0; i < len(includedlist); i++{
    for j := 0; j < len(includedlist); j++{
      if i==j{
        bndwthmat[i*len(includedlist) + j] = 1
      }
    }
  }

  for i,dp := range data_points{
    for j := 0; j< len(tmpsites); j++{
      tmpsites[j].Alpha = dp.Alpha[tmpsites[j].Aid]
    }
    // fmt.Printf("\r[TotalLike] Estimating likelihood of point %d/%d",i,len(data_points))
    all_likes[i],numno[i],totimes[i] = EstLogLikelihood(tmpsites,tmpgenes,dp.Transcript, randoms,bndwthsiz,1,bndwthmat,includedlist,rescalefactor)
  }
  loglikelihood := 0.0
  for _,lh := range all_likes{
    loglikelihood += lh
  }
  meanno,_ := stats.Mean(numno)
  //if that didn't work correctly, use jumper2 instead.
  if meanno/float64(len(randoms[0])) > 0.75{
    for i,dp := range data_points{
      for j := 0; j< len(tmpsites); j++{
        tmpsites[j].Alpha = dp.Alpha[tmpsites[j].Aid]
      }
      // fmt.Printf("\r[TotalLike] Estimating likelihood of point %d/%d",i,len(data_points))
      all_likes[i],numno[i],totimes[i] = EstLogLikelihood(tmpsites,tmpgenes,dp.Transcript, randoms,bndwthsiz,2,bndwthmat,includedlist,rescalefactor)
    }
  }
  meantim,_ := stats.Mean(totimes)
  return loglikelihood, meanno, meantim
}

//Estimate the (negative) log-likelihood of a single transctipt data point given a parameter set (including alphas)
// uses time averaging, and predrawn random numbers.
// Used in TotalLike
func EstLogLikelihood(sites []Site, genes []Gene,transcript map[string]float64, randoms [][]float64,bndwthsiz float64, whichj int, bndwthmat []float64, includedlist []int,rescalefactor float64) (float64, float64, float64) {
  //Need to get our realization into equilibrium
  //first need some initial conditions
  qdpts := 3
  num_sites := len(sites)
  num_genes := len(genes)
  var loglikest float64
  binit := make([]float64,num_sites)
  ginit := make([]float64,num_genes)
  transc_vals := make([]float64,len(includedlist))
  for i,gn := range genes{
    ginit[i] = transcript[gn.GeneID]/rescalefactor
    // ginit[i] =1
    // if ginit[i] == 0{
    //   ginit[i] = rand.Float64()
    // }
  }
  // transc_vals := make([]float64,len(includedlist))
  for i,sj := range includedlist{
    transc_vals[i] = transcript[genes[sj].GeneID]/rescalefactor//ginit[sj]
  }
  // copy(transc_vals,ginit)
  //then get a sampler going.
  path,numnoj := Jumper_use(binit,ginit,0,sites,genes,randoms,"none")
  if whichj != 1{
    path,numnoj = Jumper2_use(binit,ginit,0,sites,genes,randoms,"none")
  }
  // bndwthsiz := 2.0
  //bandwidth matrix should have det(A) = 1.


  // for i,j := range includedlist{
  //   pathg := make([]float64,len(path.T))
  //   for k :=0; k<len(path.T)-1; k++{
  //     pathg[k] = path.G[k][j]
  //   }
  //   pathavg,_ := stats.Mean(pathg)
  //   fmt.Println(transc_vals[i],pathavg,genes[j].GeneID)
  // }


  estim := 0.0

  totalT := path.T[len(path.T)-1] - path.T[1]
  all_dils := make([]float64, len(includedlist))
  for i,sj := range includedlist{
    all_dils[i] = genes[sj].Dil
  }
  // fmt.Println(all_dils)
  all_giks := make([]float64,len(includedlist))
  all_siks := make([]float64,len(includedlist))
  for j:=1; j<len(path.T)-1; j++{
    Dtk := path.T[j+1]-path.T[j]
    // copy(all_giks,path.G[j])
    for i,sj := range includedlist{
      all_giks[i] = path.G[j][sj]
    }
    // copy(all_siks,path.SG[j])
    for i,sj := range includedlist{
      all_siks[i] = path.SG[j][sj]
    }
    eval := GaussKernal(all_dils,all_giks,all_siks,transc_vals,bndwthmat,bndwthsiz)
    switch{
    case Dtk > 0:
      estim += quad.Fixed(eval, 0, Dtk, qdpts, nil, 0)
    case Dtk == 0:
      estim += 0
    case Dtk<0:
      panic("Time is moving backwards!!!")
    }
  }
  loglikest = math.Log(estim/float64(totalT))
  return -loglikest,float64(numnoj),float64(totalT)/rescalefactor
}

// Gaussian Kernal for kernal density estimate of likelihood
// Used in EstLogLikelihood, GetLiState
func GaussKernal(dis, giks, siks, data_point []float64,Aflat []float64, h float64) func(float64) float64{
  var kh float64
  return func(t float64) float64{
    model_samp := make([]float64,len(data_point))
    for i,_ := range model_samp{
      model_samp[i] = math.Exp(-dis[i]*t)*(giks[i]-siks[i]) + siks[i]
    }
    d := len(data_point)
    // fmt.Printf("d=%d, and Aflat len is %d\n",d,len(Aflat))
    Hmat := mat.NewDense(d,d,Aflat)
    msmat := mat.NewVecDense(d,model_samp)
    dpmat := mat.NewVecDense(d,data_point)

    // actual := make([]float64, d)
    // actual2 := make([]float64, d)
    xmy := mat.NewVecDense(d, nil)
    // Hxmy := mat.NewVecDense(d, actual2)

    xmy.SubVec(msmat,dpmat)
    xmy.ScaleVec(1/h,xmy)

    // Hxmy.MulVec(Hmat,xmy)
    kh = (1.0/(math.Pow(2.0*math.Pi,float64(d)/2.0)*math.Pow(h,float64(d))))*math.Exp(-0.5*mat.Inner(xmy,Hmat,xmy))
    return kh
  }
}



// FUNCTIONS TO GENERATE A REALIZATION
// ------------------------------------------------------------------------------------

//generate a realization of the model, using thinning.
// Used in GetEquilibriumSamples and main
func Jumper(b0, g0 []float64, t float64,datascl string) func([]Site,[]Gene,float64) Realization {

  g01 := make([]float64,len(g0))
  g02 := make([]float64,len(g0))
  b01 := make([]float64,len(b0))

  copy(g01, g0)
  copy(g02, g0)
  copy(b01, b0)
  state := Realization{[][]float64{b01},[][]float64{g01},[][]float64{g02},[]float64{t},[]string{}}
  return func(sites []Site, genes []Gene, endt float64) Realization {
    gids := make([]string,len(genes))
    for i,gn := range genes{
      gids[i] = gn.GeneID
    }
    state.GeneIDs = gids
    numg := len(state.G[0])
    numb := len(state.B[0])
    epign_term := make([]float64,numb)
    for i := 0; i < numb; i++{
      epign_term[i] = sites[i].Lambda*sites[i].Mu/(sites[i].Mu + math.Pow(sites[i].Alpha,sites[i].Nu))
    }
    for state.T[len(state.T) - 1] < endt {
      // PrintMemUsage()
      // runtime.GC()
      // PrintMemUsage()
      //Calculate a uniform upper bound on the sum of propensities.
      lam0 := 0.0
      cur_ind := len(state.T) - 1
      //The local steady states
      gstd := make([]float64,numg)
      for i:= 0;i<len(genes);i++{
        dot1 := 0.0
        for j,k := range genes[i].Philoc{
          dot1 += genes[i].Phival[j]*state.B[cur_ind][k]
        }
        gstd[i] = (genes[i].Gamma + dot1)/(genes[i].Dil)
      }
      copy(state.SG[cur_ind],gstd)
      //now lambda0
      gmx := make([]float64,numg)
      copy(gmx,gstd) //put max(gsteady, g(t-)) here
      for i,gs := range gmx{
        if gs < state.G[cur_ind][i]{
          gmx[i] = state.G[cur_ind][i]
        }
      }
      //now for each site add its term.
      for i:=0; i< numb; i++{
        //first compute dot(ka_i,gmx)
        dot2 := 0.0
        for j,k := range sites[i].Kaploc{
          dot2 += sites[i].Kapval[j]*gmx[k]
        }
        //then add on the thing
        lam0 += epign_term[i]*(1-state.B[cur_ind][i])*dot2 + sites[i].Lambdahat*state.B[cur_ind][i]
      }
      if lam0<0{
        panic("TIME IS MOVING BACKWARDS!!!!")
      }
      if lam0 < 0.00001{
        fmt.Println("Jumper: Next jump will not happen.")
        break
      }
      //now we have lam0. Choose the next jump time...
      // fmt.Printf("Total Propensity is %.5f \n",lam0)
      DeltaT := rand.ExpFloat64()/lam0
      //Ok we have the time to the next jump. Lets record that.
      state.T = append(state.T,state.T[cur_ind]+DeltaT)
      // fmt.Printf("Jumped to time %.5f \n",state.T[cur_ind]+DeltaT)
      //Need the value of g at (right before) the jump
      newg := make([]float64,numg)
      for i := 0; i< numg;i++{
        newg[i] = math.Exp(-(genes[i].Dil)*DeltaT)*(state.G[cur_ind][i]-state.SG[cur_ind][i]) + state.SG[cur_ind][i]
      }
      //calculate the propensities
      propens := make([]float64,2*numb)
      for i:=0; i< numb; i++{
        //first compute dot(ka_i,g)
        dot2 := 0.0
        for j,k := range sites[i].Kaploc{
          dot2 += sites[i].Kapval[j]*newg[k]
        }
        //then compute the on/off propensites
        propens[i] = epign_term[i]*(1-state.B[cur_ind][i])*dot2
        propens[i + numb] = sites[i].Lambdahat*state.B[cur_ind][i]
      }
      cumul := make([]float64,2*numb + 1)
      cumul[0] = propens[0]
      for i:= 1;i<len(cumul)-1;i++{
        cumul[i] = cumul[i-1] + propens[i]/lam0
      }
      cumul[len(cumul)-1] = 1
      uu := rand.Float64()
      choice := sort.Search(len(cumul), func(i int) bool { return cumul[i] >= uu })
      newb := make([]float64,numb)
      copy(newb,state.B[cur_ind])
      switch {
      case choice < numb://binding
        newb[choice] += 1
      case choice < 2*numb://unbinding
        newb[choice - numb] += -1
      }
      state.B = append(state.B,newb)
      state.G = append(state.G,newg)
      state.SG = append(state.SG,make([]float64,len(gstd)))
    }

  scledstate := CopyRealization(state)

  switch datascl {
  case "Log2":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log2(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log2(x)
      }
    }

  case "Log":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log(x)
      }
    }
  case "Log10":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log10(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log10(x)
      }
    }
  }

  return scledstate
  }
}

//realization that uses saved the set of random parameters, using thinning
// Used in GetEquilibriumTimeAvg,EstLogLikelihood,PrepStates
func Jumper_use(b0, g0 []float64, t0 float64,sites []Site, genes []Gene, ran_numbs [][]float64,datascl string) (Realization,int) {
  //saved set of random numbers should be array of rand.ExpFloat64 and array of rand.Float64

  numg := len(g0)
  numb := len(b0)
  num_ju := len(ran_numbs[0])
  if num_ju == 0{
    num_ju = 1
  }

  gids := make([]string,len(genes))
  for i,gn := range genes{
    gids[i] = gn.GeneID
  }
  numnojump := 0


  bvecs := make([][]float64, num_ju)
  for i := range bvecs {
    bvecs[i] = make([]float64, numb)
  }
  gvecs := make([][]float64, num_ju)
  for i := range gvecs {
    gvecs[i] = make([]float64, numg)
  }
  sgvecs := make([][]float64, num_ju)
  for i := range sgvecs {
    sgvecs[i] = make([]float64, numg)
  }



  g01 := make([]float64,len(g0))
  g02 := make([]float64,len(g0))
  b01 := make([]float64,len(b0))

  copy(g01, g0)
  copy(g02, g0)
  copy(b01, b0)

  state := Realization{bvecs,gvecs,sgvecs,make([]float64,num_ju),gids}
  state.G[0] = g01
  state.B[0] = b01
  state.SG[0] = g02
  state.T[0] = t0


  epign_term := make([]float64,numb)
  lamsum := 0.0
  for i := 0; i < numb; i++{
    epign_term[i] = sites[i].Lambda*sites[i].Mu/(sites[i].Mu + math.Pow(sites[i].Alpha,sites[i].Nu))
    lamsum += sites[i].Lambdahat
  }
  // for _,g :=range genes{
  //   fmt.Println(g.Gamma/(g.Dil))
  // }
  for ii := 0; ii < num_ju - 1; ii++{
    // fmt.Println(state.B)
    //Calculate a uniform upper bound on the sum of propensities.
    lam0 := 0.0
    //The local steady states
    gstd := make([]float64,numg)
    for i,g := range genes{
      dot1 := 0.0
      for j,k := range g.Philoc{
        dot1 += g.Phival[j]*state.B[ii][k]
      }
      gstd[i] = (g.Gamma + dot1)/(g.Dil)
    }
    //now lambda0
    gmx := make([]float64,numg)
    copy(gmx,gstd) //put max(gsteady, g(t-)) here
    for i,gs := range gmx{
      if gs < state.G[ii][i]{
        gmx[i] = state.G[ii][i]
      }
    }
    // fmt.Println(gmx)
    //now for each site add its term.
    for i,st := range sites{
      //first compute dot(ka_i,gmx)
      dot2 := 0.0
      for j,k := range st.Kaploc{
        dot2 += st.Kapval[j]*gmx[k]
      }
      //then add on the thing
      lam0 += epign_term[i]*(1-state.B[ii][i])*dot2 + sites[i].Lambdahat*state.B[ii][i]
    }
    if lam0<0{
      panic("TIME IS MOVING BACKWARDS!!!!")
    }
    if lam0 == 0{
      state.B = state.B[:ii]
      state.G = state.G[:ii]
      state.SG = state.SG[:ii]
      state.T = state.T[:ii]
      fmt.Println("Jumper_use: Next jump will not happen.")
      break
    }
    //now we have lam0. Choose the next jump time...
    ran1 := ran_numbs[0][ii]
    DeltaT := ran1/lam0
    //Ok we have the time to the next jump. Lets record that.
    state.T[ii+1] = state.T[ii]+DeltaT
    // fmt.Printf("Jumper_use: Jumping to time %0.5f\n", state.T[ii+1])
    //Need the value of g at the jump
    newg := make([]float64,numg)
    for i := 0; i< numg;i++{
      newg[i] = math.Exp(-genes[i].Dil*DeltaT)*(state.G[ii][i]-gstd[i]) + gstd[i]
    }
    // for i,_ := range newg{
    //   fmt.Printf("diff max to new: %0.5f\n",newg[i] - gmx[i])
    // }
    //calculate the propensities
    propens := make([]float64,2*numb)
    for i:=0; i< numb; i++{
      //first compute dot(ka_i,g)
      dot2 := 0.0
      for j,k := range sites[i].Kaploc{
        dot2 += sites[i].Kapval[j]*newg[k]
      }
      //then compute the on/off propensites
      propens[i] = epign_term[i]*(1-state.B[ii][i])*dot2
      propens[i + numb] = sites[i].Lambdahat*state.B[ii][i]
    }

    sumprop := 0.0
    for _,p := range propens{
      sumprop += p
    }

    // fmt.Printf("lam0: %0.5f, sum of propensites %0.5f\n", lam0,sumprop)

    cumul := make([]float64,2*numb + 1)
    cumul[0] = propens[0]/lam0
    for i:= 1;i<len(cumul)-1;i++{
      cumul[i] = cumul[i-1] + propens[i]/lam0
    }
    cumul[len(cumul)-1] = 1.0
    // fmt.Println(cumul)
    uu := ran_numbs[1][ii]
    choice := sort.Search(len(cumul), func(i int) bool { return cumul[i] >= uu })
    // fmt.Printf("num rxns %d ,choice %d\n", 2*numb ,choice)
    newb := make([]float64,numb)
    copy(newb,state.B[ii])
    switch {
    case choice < numb://binding
      newb[choice] += 1
    case choice < 2*numb://unbinding
      newb[choice - numb] += -1
    default:
      numnojump += 1
    }
    copy(state.B[ii+1],newb)
    copy(state.G[ii+1],newg)
    copy(state.SG[ii+1],gstd)
  }
  // fmt.Printf("Jumper_use: Number of non-jumps was %d\n",numnojump)
  //
  scledstate := CopyRealization(state)

  switch datascl {
  case "Log2":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log2(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log2(x)
      }
    }

  case "Log":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log(x)
      }
    }
  case "Log10":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log10(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log10(x)
      }
    }
  }

  return scledstate,numnojump
}

//generate a realization of the model that uses saved the set of random parameters, using SSA for PDMP from Crudu et al
// Faster if the bound on lambda0 is very high - leading to lots of non-jump jumps in the thinning procedure.
// Used in GetEquilibriumTimeAvg,EstLogLikelihood,PrepStates
func Jumper2_use(b0, g0 []float64, t0 float64,sites []Site, genes []Gene, ran_numbs [][]float64,datascl string) (Realization,int) {
  numg := len(g0)
  numb := len(b0)
  num_ju := len(ran_numbs[0])
  gids := make([]string,len(genes))
  for i,gn := range genes{
    gids[i] = gn.GeneID
  }

  // fmt.Println("1")
  // fmt.Println(b0)
  // fmt.Println("2")

  bvecs := make([][]float64, num_ju)
  for i := range bvecs {
    bvecs[i] = make([]float64, numb)
  }
  gvecs := make([][]float64, num_ju)
  for i := range gvecs {
    gvecs[i] = make([]float64, numg)
  }
  sgvecs := make([][]float64, num_ju)
  for i := range sgvecs {
    sgvecs[i] = make([]float64, numg)
  }


  g01 := make([]float64,len(g0))
  g02 := make([]float64,len(g0))
  b01 := make([]float64,len(b0))

  copy(g01, g0)
  copy(g02, g0)
  copy(b01, b0)

  state := Realization{bvecs,gvecs,sgvecs,make([]float64,num_ju),gids}
  state.G[0] = g01
  state.B[0] = b01
  state.SG[0] = g02
  state.T[0] = t0



  mind := genes[0].Dil

  for _,gn := range genes{
    if gn.Dil< mind{
      mind = gn.Dil
    }
  }

  if len(ran_numbs[1]) != num_ju{
    panic("[Jumper2_use] Mismatching random arrays")
  }

  state.GeneIDs = gids


  epign_term := make([]float64,numb)
  for i,site := range sites{
    epign_term[i] = site.Lambda*site.Mu/(site.Mu + math.Pow(site.Alpha,site.Nu))
  }
  for cur_ind := 0; cur_ind < num_ju-1; cur_ind++{
    // fmt.Println(state.B[cur_ind])
    // cur_ind := len(state.T) - 1 //which time point (by index) are we on
    // The local steady states
    // fmt.Println("Jumper2_use: Finding jump time")
    gstd := make([]float64,numg)
    for i:= 0;i<len(genes);i++{
      dot1 := 0.0
      for j,k := range genes[i].Philoc{
        dot1 += genes[i].Phival[j]*state.B[cur_ind][k]
      }
      gstd[i] = (genes[i].Gamma + dot1)/(genes[i].Dil)
    }
    //draw a unif[0,1]
    uu := ran_numbs[0][cur_ind]
    //we need to compute a weird little transform...
    // fmt.Println("Computing a nice little transform")
    // fmt.Println(gstd)
    stopper := FofT(state.T[cur_ind],state.G[cur_ind],state.B[cur_ind],genes,sites,gstd,epign_term)
    jumptime := state.T[cur_ind]

    qkcheck1 := stopper(math.Pow(10,5)/mind)
    if qkcheck1 < uu{
        state.B = state.B[:cur_ind]
        state.G = state.G[:cur_ind]
        state.SG = state.SG[:cur_ind]
        state.T = state.T[:cur_ind]
      fmt.Println("Jumper2_use: Next jump will not happen.")
      break
    }
    // fmt.Println("rv:",uu)
    for l := 1.0;l<8;l++{
      dt := math.Pow(10,-l)
      // fmt.Printf("Jumper2_use: random number is %0.5f\n",uu)
      for stopper(jumptime) < uu{
        // fmt.Printf("Jumper2_use: transform is %0.5f\n",stopper(jumptime))
        // fmt.Printf("Stopper is at %0.4f, needs to reach %0.4f\n",stopper(jumptime),uu)
        jumptime += dt
      }
      jumptime = jumptime-dt
    }
    if jumptime < state.T[cur_ind]{
      switch{
      case jumptime > (1.0-math.Pow(10,-8))*state.T[cur_ind]:
        jumptime = state.T[cur_ind]
      default:
        fmt.Printf("Jumper2_use: Time is not behaving. Current time is %0.6f, jumptime is somehow %0.6f\n", state.T[cur_ind], jumptime)
      }
    }
    //now the variable jumptime is when we take the jump.
    // fmt.Printf("Jumper2_use: Jumping to time %0.5f\n",jumptime)
    newg := GofT(jumptime,state.T[cur_ind],state.G[cur_ind],genes,gstd)
    //and next we have to do the normal Gillespie thing - calculate the propensities at that jump time and
    //choose from them
    propens := make([]float64,2*numb)
    for i:=0; i< numb; i++{
      //first compute dot(ka_i,g)
      dot2 := 0.0
      for j,k := range sites[i].Kaploc{
        dot2 += sites[i].Kapval[j]*newg[k]
      }
      //then compute the on/off propensites
      propens[i] = epign_term[i]*(1-state.B[cur_ind][i])*dot2
      propens[i + numb] = sites[i].Lambdahat*state.B[cur_ind][i]
    }
    lam0 := 0.0
    for _,p := range propens{
      lam0 += p
    }
    cumul := make([]float64,2*numb)
    cumul[0] = propens[0]/lam0
    for i:=1; i<len(cumul);i++{
      cumul[i] = cumul[i-1] + propens[i]/lam0
    }
    uu = ran_numbs[1][cur_ind]
    // fmt.Println("Jumper2_use: Choosing biniding/unbinding site \n")
    // fmt.Println(cumul)
    choice := sort.Search(len(cumul), func(i int) bool { return cumul[i] >= uu })
    newb := make([]float64,numb)
    copy(newb,state.B[cur_ind])
    // fmt.Printf("Jumper2_use: Choice was %d, Number of basis should be %d, Number of basis is %d\n",choice, numb, len(newb))
    switch {
    case choice < numb://binding
      newb[choice] = newb[choice] + 1
    default://unbinding
      newb[choice - numb] = newb[choice - numb] - 1
    }
    copy(state.B[cur_ind+1],newb)
    copy(state.G[cur_ind+1],newg)
    copy(state.SG[cur_ind+1],gstd)
    state.T[cur_ind+1] = jumptime
  }

  scledstate := CopyRealization(state)

  switch datascl {
  case "Log2":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log2(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log2(x)
      }
    }

  case "Log":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log(x)
      }
    }
  case "Log10":
    for i,y := range scledstate.G{
      for j,x := range y{
        scledstate.G[i][j] = math.Log10(x)
      }
    }
    for i,y := range scledstate.SG{
      for j,x := range y{
        scledstate.SG[i][j] = math.Log10(x)
      }
    }
  }


return scledstate,0
}


//  FUNCTIONS TO HELP GENERATE REALIZATIONS
// ------------------------------------------------------------------------------------

//compute -log(F) (see Crudu et al)
// Used in Jumper2_use
func FofT(tk float64, gk []float64, bk []float64, genes []Gene, sites []Site, stdy_gs []float64, epigns []float64) func(float64) float64 {
  var foft float64
  return func(t float64) float64 {
    foft = 0.0
    gintegrals := make([]float64, len(gk))
    for j,gene := range genes{
      gintegrals[j] = ((gk[j] - stdy_gs[j])/(gene.Dil))*(1-math.Exp(-gene.Dil*(t-tk))) + (t-tk)*stdy_gs[j]
    }
    kapi_d_gs := make([]float64, len(bk))
    for i,_ := range bk{
      for k,j := range sites[i].Kaploc{
        kapi_d_gs[i] += sites[i].Kapval[k]*gintegrals[j]
      }
    }

    for i,bval := range bk{
      foft += (1-bval)*epigns[i]*kapi_d_gs[i] + (t-tk)*sites[i].Lambdahat*bval
    }
    return foft
  }
}

// Compute g(t) between jumps
//used in FillIn, Jumper2_use
func GofT(t, tk float64, gk []float64, genes []Gene, stdy_gs []float64) []float64 {
  gjs := make([]float64, len(gk))
  for j,gene := range genes{
    gjs[j] = math.Exp(-(gene.Dil)*(t-tk))*(gk[j] - stdy_gs[j]) + stdy_gs[j]
  }
  return gjs
}

//get evenly spaced samples from a realization
// Used in GetEquilibriumSamples and main
func FillIn(realiz Realization, tstep float64, genes []Gene) Realization{
  fullstate := Realization{[][]float64{realiz.B[0]},[][]float64{realiz.G[0]},[][]float64{realiz.SG[0]},[]float64{realiz.T[0]}, realiz.GeneIDs}
  ct := realiz.T[0]
  for i :=0; i<len(realiz.T) - 1; i++{
    for ct <= realiz.T[i+1]{
      fullstate.B = append(fullstate.B,realiz.B[i])
      fullstate.SG = append(fullstate.SG,realiz.SG[i])
      newg := GofT(ct,realiz.T[i],realiz.G[i],genes,realiz.SG[i])
      fullstate.G = append(fullstate.G,newg)
      fullstate.T = append(fullstate.T,ct)
      ct += tstep
    }
  }
  return fullstate
}


//  FUNCTIONS TO GET EQUILIBRIUM DISTRIBUTION
// ------------------------------------------------------------------------------------
//compute the average gene transcript of the realization
// Used in GetEquilibriumSamples
func GeneAv(realiz Realization) []float64{
  roll_av := make([]float64,len(realiz.G[0]))
  for _,j := range realiz.G{
    for l,k := range j{
      roll_av[l] += k
    }
  }
  for i := 0; i < len(roll_av); i++{
    roll_av[i] = roll_av[i]/float64(len(realiz.T))
  }
  return roll_av
}

// Integrand for kernal density estimate (Gaussian)
// Used in GetEquilibriumTimeAvg
func KDEintegrand(x, di, gik, sik, h float64) func(float64) float64{
  var integrand float64
  return func(t float64) float64{
    integrand = (1/(math.Sqrt(2*math.Pi)*h))*math.Exp(-math.Pow(x-(math.Exp(-di*t)*(gik-sik) + sik),2)/(2*math.Pow(h,2)))
    return integrand
  }
}

//Returns empiricle equilibrium distrbution using time-averaging and gaussian KDE (1D)
// Used in main
func GetEquilibriumTimeAvg(sites []Site, genes []Gene, npts, qdpts ,N int,h float64,bt int, datascl string, spread float64) []EqDist{
  b0 := make([]float64,len(sites))
  g0 := make([]float64,len(genes))
  for i,_ := range(b0){
    b0[i] = float64(rand.Intn(2))
  }
  for i,gn := range(genes){
    g0[i] = (1 + gn.Gamma)/gn.Dil
  }
  jumping_ones := make([]int,0,len(genes))
  for i,gn := range genes{
    if len(gn.Philoc) > 0{
      jumping_ones = append(jumping_ones,i)
    }
  }

  //we can burn in if we really want to.
  if bt>0{

    burnrandexps := make([]float64,bt)
    burnrandunifs := make([]float64,bt)

    for i := 0; i<bt; i++{
      burnrandexps[i] = rand.ExpFloat64()
      burnrandunifs[i] = rand.Float64()
    }
    burnrandoms := [][]float64{burnrandexps,burnrandunifs}

    burnpath,_ := Jumper_use(b0, g0, 0, sites, genes, burnrandoms,"none")
    copy(b0,burnpath.B[len(burnpath.B)-1])
    copy(g0,burnpath.G[len(burnpath.G)-1])
  }

  randexps := make([]float64,N)
  randunifs := make([]float64,N)


  for i := 0; i<N; i++{
    randexps[i] = rand.ExpFloat64()
    randunifs[i] = rand.Float64()
  }

  randoms := [][]float64{randexps,randunifs}


  path,noj := Jumper_use(b0, g0, 0, sites, genes, randoms,"none")
  if float64(noj) > float64(N)*0.75{
    path,noj = Jumper2_use(b0, g0, 0, sites, genes, randoms,"none")
  }

  totalT := path.T[len(path.T)-1] - path.T[1]



  fmt.Printf("Total Realization time for estimation was %0.5f\n",totalT)

  Gtran := Transpose(path.G)
  SGtran := Transpose(path.SG)
  the_margs := make([]EqDist,len(jumping_ones))
  // for i,gn := range path.GeneIDs{
  for j,i := range jumping_ones{
    the_margs[j].GeneID = path.GeneIDs[i]

    pthmean := make([]float64,len(path.T))
    pthvar := make([]float64,len(path.T))

    pthmean[0] = Gtran[i][0]
    pthvar[0] = 0.0

    pthmean[1] = Gtran[i][1]
    pthvar[1] = 0.0



    for ti := 2; ti < len(path.T); ti ++{
      var mn float64
      var va float64
      // switch{
      // case ti<101:
      //   mn,_ = stats.Mean(Gtran[i][:ti])
      //   va,_ = stats.VarS(Gtran[i][:ti])
      // default:
      //   mn,_ = stats.Mean(Gtran[i][ti-100:ti])
      //   va,_ = stats.VarS(Gtran[i][ti-100:ti])
      // }

      mn,_ = stats.Mean(Gtran[i][:ti])
      va,_ = stats.VarS(Gtran[i][:ti])

      pthmean[ti] = mn
      pthvar[ti] = va
    }

    the_margs[j].RealMeans = pthmean
    the_margs[j].RealVars = pthvar
    the_margs[j].TimePts = path.T

    //get the domain to estimate in.
    stdev := the_margs[j].RealVars[len(path.T)-1]//stats.StdDevS(Gtran[i])
    mean :=  the_margs[j].RealMeans[len(path.T)-1]//stats.Mean(Gtran[i])
    // fmt.Println(mean,stdev)
    upper := mean + spread*stdev
    lower := mean - spread*stdev
    // h := 2.0//stdev // /float64(npts)


    switch datascl {
    case "Log2":
      the_margs[j].MeanScaled = math.Log2(mean)
      the_margs[j].StdScaled = math.Log2(stdev)
    case "Log":
      the_margs[j].MeanScaled = math.Log(mean)
      the_margs[j].StdScaled = math.Log(stdev)
    case "Log10":
      the_margs[j].MeanScaled = math.Log10(mean)
      the_margs[j].StdScaled = math.Log10(stdev)
    default:
      the_margs[j].MeanScaled = mean
      the_margs[j].StdScaled = stdev
    }



    switch datascl {
    case "Log2":
      if lower < 0{
        lower = 0.001
      }
    case "Log":
      if lower < 0.001{
        lower = 0.001
      }
    case "Log10":
      if lower < 0.001{
        lower = 0.001
      }
    default:
      if lower < 0{
        lower = 0
      }
    }



    the_margs[j].XDomain = make([]float64,npts)
    the_margs[j].DistEst = make([]float64,npts)
    delx := (upper-lower)/float64(npts)
    x := lower
    for ll :=0; ll< npts; ll+=1{
      Dx := 0.0

      switch datascl {
      case "Log2":
        the_margs[j].XDomain[ll] = math.Log2(x)
      case "Log":
        the_margs[j].XDomain[ll] = math.Log(x)
      case "Log10":
        the_margs[j].XDomain[ll] = math.Log(x)
      default:
        the_margs[j].XDomain[ll] = x
      }

      // the_margs[j].XDomain[ll] = x

      for k:=1; k<len(path.T)-1; k++{
        Dtk := path.T[k+1]-path.T[k]
        switch{
        case Dtk > 0:
          eval := KDEintegrand(x,(genes[i].Dil),Gtran[i][k],SGtran[i][k],h)
          Dx += quad.Fixed(eval, 0, Dtk, qdpts, nil, 0)
        case Dtk == 0:
          Dx += 0
        case Dtk<0:
          panic("GetEquilibriumTimeAvg: Time is moving backwards!!!")
        }
      }
      Dx = Dx/totalT
      the_margs[j].DistEst[ll] = Dx
      x += delx
    }
  }
  return the_margs
}

// Compute Transpose of a slice of floats.
// used in GetEquilibriumTimeAvg
func Transpose(slice [][]float64) [][]float64 {
    xl := len(slice[0])
    yl := len(slice)
    result := make([][]float64, xl)
    for i,_ := range result {
        result[i] = make([]float64, yl)
    }
    for i := 0; i < xl; i++ {
        for j := 0; j < yl; j++ {
            result[i][j] = slice[j][i]
        }
    }
    return result
}

//
//  FUNCTION TO GET FAKE DATA
// ------------------------------------------------------------------------------------
//create set of samples of gene expression from model - useful to make fake data
// Used in main
func GetEquilibriumSamples(b0,g0[]float64,sites []Site, genes []Gene, tol float64) func( []Site, []Gene, int) ([][]float64,[]string){
  sampler := Jumper(b0,g0,0.0,"none")
  dt := 50.0
  current_av := GeneAv(FillIn(sampler(sites,genes,10.0),0.01,genes))
  maxchng := 0.0
  for i,j := range current_av{
    if math.Abs(j-b0[i]) > maxchng{
      maxchng = math.Abs(j-b0[i])
    }
  }
  t := 20.0
  for (math.Abs(maxchng) > tol && t < 1000){
    new_av := GeneAv(FillIn(sampler(sites,genes,t),0.05,genes))
    maxchng = 0.0
    for indx,avval := range new_av{
      if math.Abs(avval-current_av[indx]) > maxchng{
        maxchng = math.Abs(avval-current_av[indx])
      }
    }
    current_av = new_av
    t+= dt
  }
  eq_after := len(FillIn(sampler(sites,genes,t),0.01,genes).T)
  //return an instance that can keep adding on samples, and even change parameters on the fly!
  // to compute an equilibrium distribution, always use nsites == sites and ngenes == genes.
  return func(nsites []Site, ngenes []Gene, num_samples int) ([][]float64,[]string){
    bigT := float64(num_samples)/100
		t+=bigT
    realization  := FillIn(sampler(nsites,ngenes,t),0.01,genes)
    sampleset := realization.G[eq_after:]
    geneids := realization.GeneIDs
    for len(sampleset) <num_samples{
      t += bigT
      sampleset = FillIn(sampler(nsites,ngenes,t),0.01,genes).G[eq_after:]
    }
    return sampleset,geneids
	}
}



func main(){


  var sitefl string
  flag.StringVar(&sitefl,"sitefl","go_in_jsons/sites.json","file name of site parameters")//path to parameters if not running from dir with them

  var genefl string
  flag.StringVar(&genefl,"genefl","go_in_jsons/genes.json","file name gene parameters")

  var choice string
  flag.StringVar(&choice,"mode","","mode: choose ParamFit, Fakeit, LogL, Eq, or none for sample path")

  var savefl string
  flag.StringVar(&savefl,"svfl","go_out_jsons/goout","name of save file")

  var datascale string
  flag.StringVar(&datascale,"DataScale","none","How the transcript data is scaled relative to raw (Log2, Log, Log10, Linear, or none). Default none. Linear is treated as none, so model will be scaled the same way. Log scaled data is converted to unscaled for all calculation, returned equilibrium and trajectories are converted back to original scale. Parameter fitting is done to unscaled data.")

  var matchedDatafl string
  flag.StringVar(&matchedDatafl,"datafl"," ","path to and name of matched data file")

  var numfake int
  flag.IntVar(&numfake,"numfake",10,"number of fake samples to produce.")

  var distres int
  flag.IntVar(&distres,"DistResolution",50,"Resolution of Equilibrium Distriubtion Estimation.")

  var qudpts int
  flag.IntVar(&qudpts,"QuadPoints",10,"Number of quadrature points for time averaging.")

  var estlength int
  flag.IntVar(&estlength,"EstimateLength",1000,"Number of jumps in log-likelihood estimator realizations.")

  var eqLength int
  flag.IntVar(&eqLength,"EqLength",1000,"Number of jumps in equilibrium distribution estimation.")

  var endtime float64
  flag.Float64Var(&endtime,"endtime",25.0,"end time of realization")

  var burnjumps int
  flag.IntVar(&burnjumps,"BurnIn",0,"Burn in for equilibrium and likelihood estimations")

  var realization_resolution float64
  flag.Float64Var(&realization_resolution,"EvenRes",0.1,"Resolution for process realizations between jumps. Set to 0 for no fill-in.")

  var whenQuit int
  flag.IntVar(&whenQuit,"StoppingCondition",10,"number of no-improvement steps before giving up")

  var toolong int
  flag.IntVar(&toolong,"maxsteps",1000,"total search steps")

  var dBand float64
  flag.Float64Var(&dBand,"DistBand",5.0,"Bandwidth for equilibrium distriubtion or likelihood estimation")

  var save_rands string
  flag.StringVar(&save_rands,"SaveRands","no","give a valid .json file name to the random parameters used in log-likelihood or parameter fitting calculation.")

  var load_rands string
  flag.StringVar(&load_rands,"LoadRands","no","give an existing .json file name to the random parameters used in log-likelihood or parameter fitting calculation.")


  var debugging_in_jupyter bool
  flag.BoolVar(&debugging_in_jupyter, "JupPrint",false,"Puts a newline at the end of each line of paramfit estimating likelihoods count, so that it runs with carriage return correctly in Jupyter for debugging.")


  var numbercores int
  flag.IntVar(&numbercores,"NumberCores",0,"Extent to parallelize making quadratic in parameter fitting (number of go routines allowed). Defauts to not parallel, use -1 for max.")

  var defaulttranscript float64
  flag.Float64Var(&defaulttranscript,"MissingTranscript",0.0,"Value to use for missing transcript in sample.")

  var defaultalpha float64
  flag.Float64Var(&defaultalpha,"MissingSite",0.0000000001,"Value to use for missing epigenetic data (alpha value in model) in sample.")


  var eq_dist_spread float64
  flag.Float64Var(&eq_dist_spread,"EqWindow",4,"Radius of equilibrium distribution to compute in units of standard deviation (i.e. number of standard deviations above & below mean to compute distribution).")

  flag.Parse()





  rand.Seed(time.Now().UnixNano())


  var thesites []Site
  var thegenes []Gene

  if (FileExists(genefl) && FileExists(sitefl)){
    thesites,thegenes = BuildModel(sitefl,genefl)
  } else {
    panic("Missing parameter files. Add path with -sitefl and -genefl flags.")
  }


  num_sites := len(thesites)
  num_genes := len(thegenes)

  fmt.Printf("Model imported with %d genes and %d binding sites\n",num_genes,num_sites)

  gene_list := make([]string,len(thegenes))
  for i,gn := range thegenes{
    gene_list[i] = gn.GeneID
  }


//site splitting means that some (real) sites can be repeated and so correspond to more than 1 site variable
  real_site := make([]string,len(thesites))
  for i,st := range thesites{
    real_site[i] = st.Aid
  }

  var data_points []MatchingData

  if FileExists(matchedDatafl){
  	datapts, _ := os.Open(matchedDatafl)
  	the_data, _ := ioutil.ReadAll(datapts)
  	json.Unmarshal(the_data,&data_points)
    fmt.Printf("Number of Data Samples found: %d \n", len(data_points))

    sites_missing := make(map[string]int)
    genes_missing := make(map[string]int)

    for _,gn := range gene_list{
      genes_missing[gn] = 0
    }

    for _,st := range real_site{
      sites_missing[st] = 0
    }


    for dpindex,dp := range data_points{

      if dp.Label == ""{
        fmt.Printf("No label for data point %d from file %s, labeling by list index\n",dpindex,matchedDatafl)
        data_points[dpindex].Label =  strconv.Itoa(dpindex)
      }


      for _, gn := range gene_list{
        if _, ok := dp.Transcript[gn]; ! ok {
          // fmt.Printf("Data is missing gene %s, assuming 0.0.\n", gn)
          genes_missing[gn] = genes_missing[gn] + 1
          dp.Transcript[gn] = defaulttranscript//we are assuming missing genes are due to low expression.
        }
      }

      track := make(map[string]bool)
      for _,st := range real_site{
        track[st] = false
      }

      for _, st := range real_site{
        if ! track[st]{
          if _, ok := dp.Alpha[st]; ! ok {
            sites_missing[st] = sites_missing[st] + 1
            dp.Alpha[st] = defaultalpha
          }
          track[st] = true
        }
      }
    }

    for _,gn := range gene_list{
      if genes_missing[gn] > 0{
        fmt.Printf("Data missing from gene %s in %d/%d samples\n",gn,genes_missing[gn],len(data_points))
      }
    }

    //Missing epigenetic data (alpha parameter) is treated as 0 by default. This has the following effect:
    //if nu >0, binding happens (and alpha>0 would reduce it)
    //if nu <0, binding WILL NOT HAPPEN.
    // Change with flag -MissingSite

    for _,st := range real_site{
      if sites_missing[st] > 0{
        fmt.Printf("Data missing from site %s in %d/%d samples\n",st,sites_missing[st],len(data_points))
      }
    }


    switch datascale{
    case "Log2":
      for _,gn :=range thegenes{
        for i,dp := range data_points{
          data_points[i].Transcript[gn.GeneID] = math.Exp2(dp.Transcript[gn.GeneID])
        }
      }
    case "Log":
      for _,gn :=range thegenes{
        for i,dp := range data_points{
          data_points[i].Transcript[gn.GeneID] = math.Exp(dp.Transcript[gn.GeneID])
        }
      }
    case "Log10":
      for _,gn :=range thegenes{
        for i,dp := range data_points{
          data_points[i].Transcript[gn.GeneID] = math.Pow(10,dp.Transcript[gn.GeneID])
        }
      }

    }




  }

  switch choice {

  //Try to optimize likelihood over parameter set
  case "ParamFit":
    fmt.Println("[ParamFit] started.")
    if FileExists(matchedDatafl){

      fmt.Println("[ParamFit] Generating random initial parameters for sites.")
      for i,_ := range thesites{
        thesites[i].Lambda = rand.Float64()
        thesites[i].Lambdahat = rand.Float64()
        thesites[i].Mu = rand.Float64()
        thesites[i].Nu = rand.Float64() - 0.5
      }
      fmt.Println("[ParamFit] Generating random initial parameters for genes.")
      for i,gn := range thegenes{
        sm := 0.0
        for _,v := range gn.Phival{
          if v < 0{
            sm = sm - v
          }
        }
        thegenes[i].Gamma = rand.Float64() + sm
        thegenes[i].Dil = rand.Float64()/100
      }


      var randoms [][]float64
      var numrands int

      switch FileExists(load_rands){
      case true:
        fmt.Println("Loading Random Numbers")
      	randop, _ := os.Open(load_rands)
      	randjs, _ := ioutil.ReadAll(randop)
      	json.Unmarshal(randjs,&randoms)
        numrands = len(randoms[0])
      case false:
        fmt.Println("Generating Random Numbers")
        numrands = estlength
        randexps := make([]float64,numrands)
        randunifs := make([]float64,numrands)
        for i := 0; i<numrands; i++{
          randexps[i] = rand.ExpFloat64()
          randunifs[i] = rand.Float64()
        }

        randoms = [][]float64{randexps,randunifs}
      }


      if strings.HasSuffix(save_rands,".json"){
        randomsave, _ := json.Marshal(randoms)
        _ = ioutil.WriteFile(save_rands, randomsave, 0644)
      }

      fmt.Printf("[ParamFit] Running MasterFit with %d randoms.\n", numrands)
      fitsites,fitgenes,fitlikelihoods := MasterFit(thesites, thegenes, data_points, 0.1,dBand, randoms,toolong,whenQuit,numbercores)

      var flnm1 string
      var flnm2 string
      var flnm3 string


      if strings.HasSuffix(savefl,".json"){
        flnm1 = savefl[:len(savefl)-5] + "_sites" + savefl[len(savefl)-5:]
        flnm2 = savefl[:len(savefl)-5] + "_genes" + savefl[len(savefl)-5:]
        flnm3 = savefl[:len(savefl)-5] + "_likelihoods" + savefl[len(savefl)-5:]
      }else{
        flnm1 = savefl + "_sites.json"
        flnm2 = savefl + "_genes.json"
        flnm3 = savefl + "_likelihoods.json"
      }

      //site parameters from fitting saved as .json
      fitsitesfl, errs := json.Marshal(fitsites)

      if errs != nil{
        fmt.Printf("Site file json encoding error: %s \n" ,errs)
      }

    	_ = ioutil.WriteFile(flnm1, fitsitesfl, 0644)

      //gene parameters from fitting saved as .json
      fitgenesfl, errg := json.Marshal(fitgenes)

      if errg != nil{
        fmt.Printf("Gene file json encoding error: %s \n" ,errg)
      }


    	_ = ioutil.WriteFile(flnm2, fitgenesfl, 0644)

      //Likelihoods computed during parameter fitting saved as .json
      fitlikesfl, _ := json.Marshal(fitlikelihoods)


    	_ = ioutil.WriteFile(flnm3, fitlikesfl, 0644)


    }else{
      fmt.Println("Data File Missing")
    }



  //Generate a bunch of fake data with fake alpha values
  case "Fakeit":
    binit := make([]float64,num_sites)
    ginit := make([]float64,num_genes)
    for i := 0; i<len(ginit); i++{
      ginit[i] = rand.Float64()
    }

    gids := make([]string,len(thegenes))
    for i,_ := range gids{
      gids[i] = thegenes[i].GeneID
    }

    alpha_id_set := make([]string,num_sites)
    for i := 0; i < num_sites; i++{
      alpha_id_set[i] = thesites[i].Aid
    }

    //numfake is a flag, default to 10
    fake_data_points := make([]MatchingData, numfake)
    for i := 0; i < numfake; i++{
      alpha_map := make(map[string]float64)
      for _,id := range alpha_id_set{
        if alpha_map[id] == 0{
          alpha_map[id] = rand.Float64()
        }
      }
      for j:= 0; j<num_sites ; j++{
        thesites[j].Alpha = alpha_map[thesites[j].Aid]
      }
      fake_data_maker := GetEquilibriumSamples(binit,ginit,thesites,thegenes,0.01)
      fdatm,geneids := fake_data_maker(thesites,thegenes,5)
      transc_map := make(map[string]float64)
      for i,gid := range geneids{
        transc_map[gid] = fdatm[len(fdatm)-1][i]
      }
      fake_data_points[i].Transcript = transc_map
      fake_data_points[i].Alpha = alpha_map
    }


    var flnm string
    if strings.HasSuffix(savefl,".json"){
      flnm = savefl
    }else{
      flnm = savefl + ".json"
    }


    fake_data_marshed,_ := json.Marshal(fake_data_points)
    _ = ioutil.WriteFile(flnm, fake_data_marshed, 0644)





  //case LogL: compute log-liklihood of a data set.
  case "LogL":
    //load matching data (data points with associated alpha values)
    if FileExists(matchedDatafl){


      var randoms [][]float64

      switch FileExists(load_rands){
      case true:
      	randop, _ := os.Open(load_rands)
      	randjs, _ := ioutil.ReadAll(randop)
      	json.Unmarshal(randjs,&randoms)

      case false:
        numrands := estlength
        randexps := make([]float64,numrands)
        randunifs := make([]float64,numrands)
        for i := 0; i<numrands; i++{
          randexps[i] = rand.ExpFloat64()
          randunifs[i] = rand.Float64()
        }

        randoms = [][]float64{randexps,randunifs}
      }



      //
      not_jumping := make([]int,0,len(thegenes))
      jumping_ones := make([]int,0,len(thegenes))
      mn_holder := make([]float64, len(data_points))
      all_means := make([]float64,len(thegenes))
      for i,gn := range thegenes{
        switch len(gn.Philoc){
        case 0:
          for j,dp := range data_points{
            mn_holder[j] = dp.Transcript[gn.GeneID]
          }
          meanof,_ := stats.Mean(mn_holder)
          all_means[i] = meanof
          not_jumping = append(not_jumping,i)
        default:
          for j,dp := range data_points{
            mn_holder[j] = dp.Transcript[gn.GeneID]
          }
          meanof,_ := stats.Mean(mn_holder)
          all_means[i] = meanof
          jumping_ones = append(jumping_ones,i)
        }
      }

      final_mean,_ := stats.Mean(all_means)



      if strings.HasSuffix(save_rands,".json"){
        randomsave, _ := json.Marshal(randoms)
        _ = ioutil.WriteFile(save_rands, randomsave, 0644)
      }




      for i,gn := range thegenes{
        thegenes[i].Dil = gn.Dil * final_mean
      }







      loglikelihood,nojumps,avgtimes := TotalLike(thesites,thegenes,data_points,randoms,dBand,jumping_ones,final_mean)

      fmt.Printf("Log-Likelihood of was %0.5f, average simulation did not jump on %0.5f/%d attempts, average simulation time length was %0.5f\n",loglikelihood,nojumps,estlength, avgtimes)
    }else{
      fmt.Println("Data File Missing")
    }





  //case Eq: get equilibrium distribution samples
  case "Eq":
    if FileExists(matchedDatafl){


      for _,dp := range data_points{
        for j,_ := range thesites{
          thesites[j].Alpha = dp.Alpha[thesites[j].Aid]
        }
        eq_samples := GetEquilibriumTimeAvg(thesites,thegenes,distres,qudpts,eqLength,dBand,burnjumps,datascale,eq_dist_spread)



        for j,gen := range eq_samples{
          for k,val := range gen.DistEst{
            if math.IsInf(val,0){
              eq_samples[j].DistEst[k] = -1.0
              fmt.Printf("%s has Inf value\n",gen.GeneID)
            }
            if math.IsNaN(val){
              eq_samples[j].DistEst[k] = -2.0
              fmt.Printf("%s has NaN value\n",gen.GeneID)
            }
          }
          for k,val := range gen.RealMeans{
            if math.IsInf(val,0){
              eq_samples[j].RealMeans[k] = -1.0
              fmt.Printf("%s has Inf value\n",gen.GeneID)
            }
            if math.IsNaN(val){
              eq_samples[j].DistEst[k] = -2.0
              fmt.Printf("%s has NaN value\n",gen.GeneID)
            }
          }
          for k,val := range gen.RealVars{
            if math.IsInf(val,0){
              eq_samples[j].RealVars[k] = -1.0
              fmt.Printf("%s has Inf value\n",gen.GeneID)
            }
            if math.IsNaN(val){
              eq_samples[j].DistEst[k] = -2.0
              fmt.Printf("%s has NaN value\n",gen.GeneID)
            }
          }
        }




        eq_file,err := json.Marshal(eq_samples)


        if err != nil{
          fmt.Printf("Sample %s json encoding error: %s \n", dp.Label ,err)
        }


        var flnm string
        if strings.HasSuffix(savefl,".json"){
          flnm = savefl[:len(savefl)-5] + "_"+dp.Label + ".json"
        }else{
          flnm = savefl + "_"+dp.Label+".json"
        }


        _ = ioutil.WriteFile(flnm, eq_file, 0644)
      }
    }else{
      eq_samples := GetEquilibriumTimeAvg(thesites,thegenes,distres,qudpts,eqLength,dBand,burnjumps,datascale,eq_dist_spread)


      for j,gen := range eq_samples{
        for i,val := range gen.DistEst{
          if math.IsInf(val,0){
            eq_samples[j].DistEst[i] = -1.0
          }
        }
      }

      for j,gen := range eq_samples{
        for i,val := range gen.RealMeans{
          if math.IsInf(val,0){
            eq_samples[j].RealMeans[i] = -1.0
          }
        }
      }

      for j,gen := range eq_samples{
        for i,val := range gen.RealVars{
          if math.IsInf(val,0){
            eq_samples[j].RealVars[i] = -1.0
          }
        }
      }

      eq_file,err := json.Marshal(eq_samples)

      if err != nil{
        fmt.Printf("json encoding error: %s \n" ,err)
      }

      var flnm string
      if strings.HasSuffix(savefl,".json"){
        flnm = savefl
      }else{
        flnm = savefl + ".json"
      }


      _ = ioutil.WriteFile(flnm, eq_file, 0644)
    }







  //default case: get a realization of the process
  default:


      fmt.Println("Making trajectories")


      binit := make([]float64,len(thesites))
      for j,_ := range(binit){
        binit[j] = float64(rand.Intn(2))
      }

      ginit := make([]float64,len(thegenes))

      // var preinit float64

      for i,gn := range(thegenes){
        ginit[i] = ((gn.Gamma + 1)/(gn.Dil))*(0.75 + 0.5*rand.Float64())
      }

      fmt.Println(ginit)



      if FileExists(matchedDatafl){

        for _,dp := range data_points{

          for j,_ := range thesites{
            thesites[j].Alpha = dp.Alpha[thesites[j].Aid]
          }

          runnit := Jumper(binit,ginit,0.0,datascale)
          resul := runnit(thesites,thegenes,endtime)


          for k,ar := range resul.G{
            for j,val := range ar{
              if math.IsInf(val,0){
                resul.G[k][j] = -1.0
              }
            }
          }


          for k,ar := range resul.SG{
            for j,val := range ar{
              if math.IsInf(val,0){
                resul.SG[k][j] = -1.0
              }
            }
          }


          var flnm1 string
          var flnm2 string
          if strings.HasSuffix(savefl,".json"){
            flnm1 = savefl[:len(savefl)-5] + "_jumps_"+ dp.Label + savefl[len(savefl)-5:]
            flnm2 = savefl[:len(savefl)-5] + "_even_" + dp.Label + savefl[len(savefl)-5:]
          }else{
            flnm1 = savefl + "_jumps_"+ dp.Label+ ".json"
            flnm2 = savefl + "_even_"+dp.Label +".json"
          }

          realizfile, err := json.Marshal(resul)

          if err != nil{
            fmt.Printf("Sample %s json encoding error: %s \n", dp.Label ,err)
          }


        	_ = ioutil.WriteFile(flnm1, realizfile, 0644)


          if realization_resolution > 0{
            fullres := FillIn(resul,realization_resolution,thegenes)
            fulrealizfile, errf := json.Marshal(fullres)


            if errf != nil{
              fmt.Println("json encoding error (full realization):", errf)
            }


          	_ = ioutil.WriteFile(flnm2, fulrealizfile, 0644)

          }


        }
      }else{


        ginit := make([]float64,len(thegenes))

        for i,gn := range(thegenes){
          thegenes[i].Dil = gn.Dil
          ginit[i] = (gn.Gamma/(gn.Dil))*(0.75 + 0.5*rand.Float64())
        }

        runnit := Jumper(binit,ginit,0.0,"none")
        resul := runnit(thesites,thegenes,endtime)

        for i,ar := range resul.G{
          for j,val := range ar{
            if math.IsInf(val,0){
              resul.G[i][j] = -1.0
            }
          }
        }


        for i,ar := range resul.SG{
          for j,val := range ar{
            if math.IsInf(val,0){
              resul.SG[i][j] = -1.0
            }
          }
        }



        var flnm1 string
        var flnm2 string
        if strings.HasSuffix(savefl,".json"){
          flnm1 = savefl[:len(savefl)-5] + "_jumps" + savefl[len(savefl)-5:]
          flnm2 = savefl[:len(savefl)-5] + "_even" + savefl[len(savefl)-5:]
        }else{
          flnm1 = savefl + "_jumps.json"
          flnm2 = savefl + "_even.json"
        }

        realizfile, err := json.Marshal(resul)
        if err != nil{
          fmt.Println("json encoding error:", err)
        }

      	_ = ioutil.WriteFile(flnm1, realizfile, 0644)

        if realization_resolution > 0{
          fullres := FillIn(resul,realization_resolution,thegenes)

          fulrealizfile, errf := json.Marshal(fullres)


          if errf != nil{
            fmt.Println("json encoding error (full realization):", errf)
          }


        	_ = ioutil.WriteFile(flnm2, fulrealizfile, 0644)
        }


      }
   }

}
