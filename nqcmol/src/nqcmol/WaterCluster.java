/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

/**
 *
 * @author nqc
 */
public class WaterCluster extends Cluster {

	final private int cSigMax=5;

	protected int[] sigO=new int [cSigMax];

	/**
	 * Get the value of sigO
	 *
	 * @return the value of sigO
	 */
	public int[] getSigO() {
		return sigO;
	}

	/**
	 * Get the value of sigO at specified index
	 *
	 * @param index
	 * @return the value of sigO at specified index
	 */
	public int getSigO(int index) {
		return this.sigO[index];
	}


	boolean CalcSigO(int verbose){
		int j,count;
		double dist2;
		for(int i=0;i<cSigMax;i++) sigO[i]=0;
		for(int i=0;i<nAtoms;i++)	if(Nz[i]==8)
			for(int j=0;j<nAtoms;j++)	if(Nz[j]==8){
				if(IsConnected(i,j))
				count=nOO[i];
				if(count>=cSigMax){
					cout<<" A strange structure appears ! Cannot identify ! "<<endl;
					return -1;
				}else sigO[count]++;
			}
			return true;
	};

	int CalcOH(int verbose){
		int i,j,answer=0;
		//cout<<"OH ="<<bondOH.t1<<" "<<bondOH.t2<<" "<<bondOH.dmin<<" "<<bondOH.dmax<<endl;
		//cout<<"HO ="<<bondHO.t1<<" "<<bondHO.t2<<" "<<bondHO.dmin<<" "<<bondHO.dmax<<endl;
		for(i=0;i<nAtoms;i++) nOH[i]=nOO[i]=0;
		for(i=0;i<nAtoms;i++){
			//cout<<"xHO ="<<Nz[i];
			for(j=0;j<nAtoms;j++) if(i!=j){
				//cout<<" "<<Nz[j]<<" "<<distance(i,j)<<endl;

				if(bBond2(i,j,bondOH)){ //get all neighbor H atoms  -- for O
					if(nOH[i]<4){	lOH[i][nOH[i]]=j;		nOH[i]++;
					}else{ answer=-1; goto verbose_calcOH;};
				}

				if(bBond2(i,j,bondHO)){ //get all neighbor O atoms -- for H, the covalent bond put first
					if(nOH[i]<2){
						lOH[i][nOH[i]]=j;		nOH[i]++;
						if((nOH[i]==2)&&(distance(i,lOH[i][0])>distance(i,lOH[i][1])))
							swap(lOH[i][0],lOH[i][1]);
					}else{
						answer=-2;
						goto verbose_calcOH;
					};
				}
			}
		}

	//get all neighbor H atoms  -- for O
	for(i=0;i<nAtoms;i++) if((Nz[i]==HYDROGEN)&&(nOH[i]==2)){
		int l1=lOH[i][0],l2=lOH[i][1];
		if(nOO[l1]<4){ lOO[l1][nOO[l1]]=l2;nOO[l1]++;}else{ answer=-3; goto verbose_calcOH;};

		if(nOO[l2]<4){ lOO[l2][nOO[l2]]=l1;nOO[l2]++; }else{ answer=-3; goto verbose_calcOH;};
	}

	verbose_calcOH:
	if(verbose){
		cout<<"Answer = "<<answer<<endl;
		for(i=0;i<nAtoms;i++){
			cout<<" Atom "<<i+1<<" has "<<nOH[i]<<" O-H-bonds : ";
			for(j=0;j<nOH[i];j++)
				cout<<lOH[i][j]+1<<" ";
			cout<<endl;
		}

		for(i=0;i<nAtoms;i++) if(Nz[i]==OXYGEN){
			cout<<" Oxygen "<<i+1<<" has "<<nOO[i]<<" O-O-bonds : ";
			for(j=0;j<nOO[i];j++)
				cout<<lOO[i][j]+1<<" ";
			cout<<endl;
		}
	}
	return answer;
}



	int getNumberOfSmallestRing(){
		CalcSigO();
		int d=0;
		for(int i=0;i<cSigMax;i++){
			d+=sigO[i]*i;
		}
		return d/2 - nType[8] +1;
	}

}
