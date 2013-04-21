/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package joey.Structure;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import org.biojava.bio.structure.Atom;
import org.biojava.bio.structure.Structure;
import org.biojava.bio.structure.StructureTools;
import org.biojava.bio.structure.io.PDBFileReader;
import org.omg.PortableInterceptor.SYSTEM_EXCEPTION;

/**
 *
 * @author zhongxf
 */
public class ParseStructure 
{
    double[][] rmsdMatrix; // hold the rmsd of structures
    String pwdString = System.getProperty("user.dir"); //hold current work directory

    public void ParsStruct(Structure[] structArray) 
    {
        int index =0;
        calculateRmsd(index,structArray);
        rmsdMatrix = calculateRmsd(structArray);
    }
    public double[][] getRmsdMatrix()
    {
        return rmsdMatrix;
    }

    private void calculateRmsd(int index, Structure[] structArray) 
    {
        
    }
    private double[][] calculateRmsd(Structure[] structArray) 
    {
        double[][] RmsdArray = new double[structArray.length][structArray.length];      
        Atom[][]  coordsArray = new Atom[structArray.length][];
        //initial
        for (int i =0; i < structArray.length; i++)
        {
            coordsArray[i] = StructureTools.getAtomCAArray(structArray[i]);
        }
        
        int coordsArrayLength = coordsArray.length;
        
        for(int j=0; j< coordsArrayLength;j++)
        {
            int atomLen = coordsArray[j].length;           
            double sumTmp =0;           
            for (int k = 0; k < coordsArrayLength; k++)
            { 
                if (k == j) 
                {
                    RmsdArray[j][k] = 0;
                } 
                else 
                {
                    if (k < j) 
                    {
                        RmsdArray[j][k] = RmsdArray[k][j];                       
                    } 
                    else 
                    {
                        for (int m = 0; m < atomLen; m++) 
                        {
                            try 
                            {
                                double[] coords1 = null;
                                double[] coords2 = null;;
                                coords1 = coordsArray[j][m].getCoords();
                                coords2 = coordsArray[k][m].getCoords();                               
                                double atomRmsd = (Math.pow((coords1[0] - coords2[0]), 2)
                                        + Math.pow((coords1[1] - coords2[1]), 2)
                                        + Math.pow((coords1[2] - coords2[2]), 2));

                                sumTmp += atomRmsd;
                            } 
                            catch (Exception ex) 
                            {
                                ex.printStackTrace();
                            }                           
                        }
                         RmsdArray[j][k] = Math.sqrt( sumTmp / atomLen);
                    }
                }               
            }
        }  
        return RmsdArray;
    }

    public Structure[] parsFileArray(File[] files) throws FileNotFoundException, IOException
    {
        ArrayList struArrayList = new ArrayList();
        for (int i =0; i<files.length; i++)
        {
            Structure[] tmpStructures =parsFile(files[i]);
            for (int j=0; j<tmpStructures.length; j++)
            {
                struArrayList.add(tmpStructures[j]);
            }
        }
        //ArrayList convert to Array
        Structure[] struArray = new Structure[struArrayList.size()];
        for (int k=0; k<struArrayList.size(); k++)
        {
            struArray[k] = (Structure)struArrayList.get(k);
        }
        return struArray;
    }
    public Structure[] parsFile(File file) throws FileNotFoundException, IOException 
    {
        String readLineString;
        FileReader fileReader = new FileReader(file);
        BufferedReader buffReader = new BufferedReader(fileReader);
        ArrayList struStrList = new ArrayList();
        readLineString = buffReader.readLine();
        String tmpString = "";
        while (readLineString != null)
        {           
            if (readLineString.startsWith("ATOM"))
            {
                tmpString +=(readLineString + "\n");
            }
            else 
            {
                if (!tmpString.equals(""))
                {
                    struStrList.add(tmpString);
                }
                tmpString = "";
            }
            readLineString = buffReader.readLine();
        }
        File[] struFiles = writeToFiles(struStrList);
        Structure[] structArray = extractFromFiles(struFiles);
        return structArray;
    }

    private File[] writeToFiles(ArrayList struStrList) throws IOException 
    {
        File[] struFiles = new File[struStrList.size()];
        for (int i=0; i < struStrList.size(); i++)
        {
//            File tmpFile = new File(pwdString+"/test" + (i + 1) + ".pdb");
            File tmpFile = new File("/home/zhongxf/atest/test/data0407/test" + (i + 1) + ".pdb");
            if (tmpFile.exists())
            {
                System.out.println("File exists");
            }
            else 
            {
                if (tmpFile.createNewFile())
                {
                    BufferedWriter bufferedWriter = new BufferedWriter(new FileWriter(tmpFile));
                    bufferedWriter.write((String) struStrList.get(i));
                    bufferedWriter.flush();
                    bufferedWriter.close();
                }
            }
            struFiles[i] = tmpFile;
        }
        return struFiles;
    }

    private Structure[] extractFromFiles(File[] struFiles) throws IOException 
    {
        Structure[] structures = new Structure[struFiles.length];
        for (int i =0; i < struFiles.length; i++)
        {
            PDBFileReader pdbRead = new PDBFileReader();
            File tmpFile = struFiles[i];
            System.out.println(tmpFile.toString());
            Structure tmpStru = pdbRead.getStructure(tmpFile);
            structures[i] =  tmpStru;
        }
        return structures;
    }

}
