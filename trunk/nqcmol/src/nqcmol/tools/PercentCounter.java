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
public class PercentCounter {
    HashMap<String,Integer> data=new HashMap<String,Integer>();

    void increment(String k){
         if( data.containsKey(k) ){
            data.put(k, data.get(k)+1);
            //cout<<" k ++"<<k<<" "<<ValSpan(k)<<" "<<bin[k]<<" "<<bin.max_size()<<" "<<bin.size()<<endl;
        }else	data.put(k,1);
    }

    void getValue(String k){

    }
}
