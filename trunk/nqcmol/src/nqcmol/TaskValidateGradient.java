/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;


/**
 *
 * @author nqc
 */
public class TaskValidateGradient extends TaskCalculate {
	@Override
	public String getName(){
		return "ValidateGradient";
	}

	@Override
	protected void Process() {		
//			xmllog.writeEntity("Note");
//			xmllog.writeText(" Time is measured in seconds");
//			xmllog.endEntity();
			
			int i = 0;
			while (mol.Read(fileIn,sFormatIn)) {
				//mol.Write(System.out,"xyz");
				pot.setCluster(mol);
				System.out.print(pot.ValidateGradient(1.0E-4));

				i++;
			}
			fileIn.close();				
	}

	
}
