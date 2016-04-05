import java.util.*;
import java.io.*;

/**
*proLigContacts - lists resiudes in contact with biologically relevant ligands in a PDB file
*/
public class proLigContacts
{
	Hashtable bindingResiduesHash = new Hashtable();
	Hashtable bindingResiduesTypeHash = new Hashtable();
	Hashtable ligandTypeHash = new Hashtable();
	Hashtable resultshash = new Hashtable();
	
	public proLigContacts( String pdbfile )
	{
		try
		{
			//parse Mark's list and add to relevantLigandSet
			//System.out.println(pdbfile);
// 			/*/*HashSet relevantLigandSet = new HashSet();
// 			InputStream is1 = proLigContacts.class.getResourceAsStream("/uniprotLigandListExtended_classified.txt");
// 			BufferedReader in1 = new BufferedReader(new InputStreamReader(is1));
// 			String ligline = in1.readLine();
// 			while (null != (ligline = in1.readLine()))
// 			{
// 				if( !ligline.startsWith( "#" ) && !ligline.startsWith( "//" ) )
// 				{
// 					if( ligline.length() > 0 )
// 					{
// 						//System.out.println( ligline );
// 						String ligandabbr =  ligline.substring( 0, ligline.indexOf("=") );
// 						//System.out.println( ligandabbr );
// 						relevantLigandSet.add( ligandabbr );
// 						
// 					}
// 				}
// 			}
// 			in1.close();*/*/
			//System.out.println( relevantLigandSet );
			
			//System.out.println("HERE");

			//check HETNAM records for those that agree with fireDB output
			//get coords for records that agree
			if( ( new File( pdbfile ) ).exists() )
			{
				
				//System.out.println( "Finding coords for biologically relevant ligands in " + pdbfile + "..."); TEMP
				BufferedReader in = new BufferedReader( new FileReader( pdbfile ) );
				String line = in.readLine();
				StringBuffer coordsbuf = new StringBuffer();
				boolean relevant = false;
				do
				{
					if( line.startsWith( "ATOM" ) || line.startsWith( "HETATM" ) )
					{
						//System.out.println( line );
						
						//if( line.charAt( 21 ) ==  chain )
						//{
							if( line.startsWith( "ATOM" ) )
							{
							      coordsbuf.append(line+"\n");
							}
							
							//get type of HETATM only output coords of relevant ligands
							if( line.startsWith( "HETATM" ) )
							{
								//get type of HETATM
								String het_type = line.substring( 17, 20 ).trim();
								//System.out.println(het_type);
								
								//if( relevantLigandSet.contains( het_type ) )
								//{
									coordsbuf.append(line+"\n");
									relevant = true;
								//}
							}
						//}	
					}
					
					line = in.readLine();
				}
				while( line != null && !line.startsWith( "ENDMDL" ) && !line.startsWith( "MODEL        2" ) );
				in.close();
				
				//System.out.println("NOW HERE");
				//System.out.println(coordsbuf);
				
				if( relevant )
				{
					//output coords to a new PDB file
					DataOutputStream out = new DataOutputStream( new FileOutputStream( pdbfile + "_lig_temp.pdb" ) );
					out.writeBytes( coordsbuf.toString() );
					//System.out.println( coordsbuf.toString() );
					out.close();
					//System.out.println("TEST 1234");
					
					//make a Hashtable containing Van der Waals radii info
					Hashtable vdwHash = new Hashtable();
					InputStream is = proLigContacts.class.getResourceAsStream("/vanderwaalsradii.dat");
					in = new BufferedReader(new InputStreamReader(is));
					line = in.readLine();
					//System.out.println(line);
					
					while (null != (line = in.readLine()))
					{
						if( line.length() > 0 )
						{
							StringTokenizer linetoks = new StringTokenizer( line, "," );
							String element = linetoks.nextToken().toUpperCase();
							Float radius = ( new Float( linetoks.nextToken() ) );
							vdwHash.put(element, radius);
						}
					}
					in.close();
					is.close();
					
					//store coordinates of ATOMS and HETATMS
					in = new BufferedReader( new FileReader( pdbfile + "_lig_temp.pdb" ) );
					line = in.readLine();
					
					//System.out.println("TEST XYZ");
					
					DataOutputStream out3 = new DataOutputStream( new FileOutputStream( pdbfile + "_hetatm.out" ) );
					
					StringBuffer ATOMcoordsbuf = new StringBuffer();
					StringBuffer HETATMcoordsbuf = new StringBuffer();
					do
					{
						if( line.startsWith( "ATOM" ) )
							ATOMcoordsbuf.append(line+"\n");

						if( line.startsWith( "HETATM" ) )
							HETATMcoordsbuf.append(line+"\n");
							out3.writeBytes(line+"\n");
							
						
						line = in.readLine();
					}
					while( line != null );
					in.close();

					//System.out.println("TEST ABC");
					
					//calculate distances between all HETATMs and residue ATOMs (compare HETATM coords with ATOM coords)
					String ATOMcoordsstr = ATOMcoordsbuf.toString();
					String HETATMcoordsstr = HETATMcoordsbuf.toString();
					StringTokenizer HETATMtoks = new StringTokenizer( HETATMcoordsstr, "\n" );
					
					while(HETATMtoks.hasMoreTokens())
					{
						String line1 = HETATMtoks.nextToken();
						//System.out.println(line1);
						String atominf = line1.substring( 0, 26 );
						//System.out.println(atominf);
						String atomtype = atominf.substring(13, 14).trim();
						//String atomtype = line1.substring(76).trim();
						
						//System.out.println( " atomtype " +atomtype );
						
						//handle elements in new version of pymol e.g. O1- instead of just O
						StringBuffer atomtypebuf = new StringBuffer();
						for( int e = 0; e < atomtype.length(); e++ )
						{
							char atomchar =  atomtype.charAt( e );
							if( Character.isLetter( atomchar ) )
							{
								atomtypebuf.append( atomchar );
								//System.out.println(" atomtypebuf " + atomtypebuf);
							}
						}
						atomtype = atomtypebuf.toString();
// 				
						String ligtype = atominf.substring(17, 20).trim();
						//String ligtype = atominf.substring(76).trim();
						//System.out.println(" ligtype "+ligtype);
						int atomnum = (new Integer(atominf.substring(7, 12).trim()));
						//System.out.println("atomnum " + atomnum);
						int lignum = (new Integer(atominf.substring(22, 26).trim()));
						//System.out.println("lignum " + lignum);
						String xstr = line1.substring(30,38).trim();
						String ystr = line1.substring(38,46).trim();
						String zstr = line1.substring(47,54).trim();
						//System.out.println("xstr " + xstr);
						//System.out.println("ystr " + ystr);
						//System.out.println("zstr " + zstr);
						
						
						if( !xstr.equals( "nan" ) && !ystr.equals( "nan" ) && !zstr.equals( "nan" ) )
						{
							float x1 = (new Float(xstr)).floatValue();
							float y1 = (new Float(ystr)).floatValue();
							float z1 = (new Float(zstr)).floatValue();
							
							StringTokenizer ATOMtoks = new StringTokenizer( ATOMcoordsstr, "\n" );
							
							while(ATOMtoks.hasMoreTokens())
							{
								String line2 = ATOMtoks.nextToken();
								//System.out.println("line2 " + line2 );
								String atominf2 = line2.substring( 0, 26 );
								String atomtype2 = atominf2.substring(13, 14).trim();
								//String atomtype2 = line2.substring(76).trim();
								
								//System.out.println( "atomtype2 " + atomtype2 );
								
								//handle elements in new version of pymol e.g. O1- instead of just O
								StringBuffer atomtype2buf = new StringBuffer();
								for( int e = 0; e < atomtype2.length(); e++ )
								{
									char atomchar =  atomtype2.charAt( e );
									if( Character.isLetter( atomchar ) )
									{
										atomtype2buf.append( atomchar );
										//System.out.println("atomtype2buf " +atomtype2buf);
									}
								}
								atomtype2 = atomtype2buf.toString();
								
								//String restype2 = atominf2.substring(76).trim();
								String restype2 = atominf2.substring(17, 20).trim();
								//System.out.println("restype2 "+restype2);
								int atomnum2 = (new Integer(atominf2.substring(7, 12).trim()));
								//System.out.println("atomnum2 "+atomnum2);
								int resnum2 = (new Integer(atominf2.substring(22, 26).trim()));
								//System.out.println("resnum2 "+resnum2);
								String xstr2 = line2.substring(30,38).trim();
								//System.out.println("xstr2 "+xstr2);
								String ystr2 = line2.substring(38,46).trim();
								//System.out.println("ystr2 "+ystr2);
								String zstr2 = line2.substring(47,54).trim();
								//System.out.println("zstr2 "+zstr2);
								
								if( !xstr2.equals( "nan" ) && !ystr2.equals( "nan" ) && !zstr2.equals( "nan" ) )
								{
									float x2 = (new Float(xstr2)).floatValue();
									float y2 = (new Float(ystr2)).floatValue();
									float z2 = (new Float(zstr2)).floatValue();
									
									//calc Euclidean distance between HETATM and ATOM
									double distance = Math.sqrt( ((x1-x2)*(x1-x2)) + ((y1-y2)*(y1-y2)) + ((z1-z2)*(z1-z2)) );
									//System.out.println(" distance " + distance);
									
									//get vdw radii
									float distHETATM = ((Float)vdwHash.get( atomtype )).floatValue();
									//System.out.println("distHETATM  " + distHETATM);
									float distATOM = ((Float)vdwHash.get( atomtype2 )).floatValue();
									//System.out.println("distATOM " + distATOM);
									
									if( distance <= distHETATM + distATOM + 0.5 )
									{
										//System.out.println( ligtype + " - " + atomtype + " (" + distHETATM + ") in contact with " + atomtype2 + " (" + distATOM + "): " + distance + " Angstroms\n" + atominf + "\n" + atominf2 );
										//System.out.println(" distance " + distance);
										
										//store info on binding residues
										if( !bindingResiduesHash.containsKey( new Integer( lignum ) ) )
										{
											TreeSet residueSet = new TreeSet();
											//String resi = resnum2+" ";//;+restype2;
											//residueSet.add( resi );
											residueSet.add( resnum2 );
											bindingResiduesHash.put( new Integer( lignum ), residueSet );
/*											
											TreeSet residueType = new TreeSet();
											residueType.add( restype2 );
											bindingResiduesTypeHash.put( new Integer( lignum ), residueType );
											*/
											//System.out.println("TEST ABC" + resnum2 + " " + restype2);
										}
										
										else if( bindingResiduesHash.containsKey( new Integer( lignum ) ) )
										{
											//only count ligand-residue interaction once for same ligand
											TreeSet residueSet = (TreeSet)bindingResiduesHash.get( new Integer( lignum ) );
											Integer resi = resnum2;
											//String resi = resnum2+" ";//+restype2;
											residueSet.add( resi );
											residueSet.add( resnum2 );
											bindingResiduesHash.put( new Integer( lignum ), residueSet );
											
// 											TreeSet residueType = (TreeSet)bindingResiduesTypeHash.get( new Integer( lignum ) );
// 											residueType.add( restype2 );
// 											bindingResiduesTypeHash.put( new Integer( lignum ), residueType );
											
											//System.out.println(resnum2 + " " + restype2);
											//System.out.println( "test 123 " + resnum2 + " " + restype2);
										}

										//store info on ligand types
										ligandTypeHash.put( new Integer( lignum ), ligtype );
										
									}
								}
								
							}
						}
					}

					//System.out.println( "TEST END" );
					//System.out.println( bindingResiduesHash );
					System.out.println( "\nType\tBinding residues"); 
					//String pdbid =  pdbfile.substring( 0, pdbfile.indexOf("_") );
					//String pdbid = pdbfile;
					DataOutputStream out2 = new DataOutputStream( new FileOutputStream(  pdbfile+ "_bs.out" ) );
					
					for( Enumeration enumer = bindingResiduesHash.keys(); enumer.hasMoreElements(); )
					{
						Integer lignum = (Integer) enumer.nextElement();
						TreeSet residueSet = (TreeSet)bindingResiduesHash.get( new Integer( lignum ) );
						//System.out.println(residueSet);
						String ligtype = (String)ligandTypeHash.get( new Integer( lignum ) );
						System.out.println( ligtype + "\t" + residueSet );  //TEMP

						out2.writeBytes( ligtype + "\t" + residueSet  );
						resultshash.put( residueSet.toString(), ligtype ); 

					}
					out2.close();
					
					//System.out.println(" plc resultshash " + resultshash  );
					
					
				}

				else
				{
					System.out.println( "No biologically relevant ligands found." );

				}
			}
		}

		catch( Exception e )
		{
			System.err.println( "Error executing proLigContacts!\n" +e );
		}
		
	}
	
	public Hashtable getresultshash()
	{
		return resultshash;
	}


	public static void main( String args[])
	{
		try
		{
			proLigContacts plc = new proLigContacts( args[0] );
			Hashtable resultshash = (Hashtable)plc.getresultshash();
			//System.out.println("plc test " + resultshash);

		}
		catch( Exception e )
		{
			System.err.println( e );
		}
	}
}
