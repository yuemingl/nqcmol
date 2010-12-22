/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.tools;

import java.util.HashMap;


/**
 *
 * @author chinh
 */
public class Histogram extends org.apache.commons.math.stat.descriptive.DescriptiveStatistics{

    @Override
    public void addValue(double v) {
        int k=(int)((v-refValue)/binSize);
        //cout<<" k ++"<<k<<" "<<ValSpan(k)<<endl;
        if( bin.containsKey(k) ){
            bin.put(k, bin.get(k)+1);
            //cout<<" k ++"<<k<<" "<<ValSpan(k)<<" "<<bin[k]<<" "<<bin.max_size()<<" "<<bin.size()<<endl;
        }else	bin.put(k,1);

        super.addValue(v);
    }

    protected int refValue = 0;

    /**
     * Get the value of refValue
     *
     * @return the value of refValue
     */
    public int getRefValue() {
        return refValue;
    }

    /**
     * Set the value of refValue
     *
     * @param refValue new value of refValue
     */
    public void setRefValue(int refValue) {
        this.refValue = refValue;
    }


    protected int binSize = 1;

    public int getBinSize() {
        return binSize;
    }

    public void setBinSize(int binSize) {
        this.binSize = binSize;
    }

    HashMap<Integer,Integer> bin=new HashMap<Integer,Integer>();
}
