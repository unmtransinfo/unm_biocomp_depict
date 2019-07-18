package edu.unm.health.biocomp.depict;

import java.io.*;
import java.text.*;
import java.util.*;
import java.util.regex.*;
import java.awt.Color;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.util.MolHandler;
import chemaxon.license.LicenseManager;


/**	Utility/test app: Input one molecule, generates PNG or JPEG image.
	<br>
	@author Jeremy J Yang
*/
public class mol2img_app
{
  private static int verbose=0;
  private static String ifile=null;
  private static String ofile=null;
  private static Integer width=300;
  private static Integer height=240;
  private static String imgfmt="PNG";
  private static Integer maxscale=28;
  private static Boolean kekule=false;
  private static Boolean arom_gen=false;
  private static Boolean arom_bas=false;
  private static Boolean showh=false;
  private static Boolean transparent=false;
  private static String smarts=null;

  private static void Help(String msg)
  {
    System.err.println(msg+"\n"
      +"mol2img - depict utility\n"
      +"\n"
      +"usage: mol2img [options]\n"
      +"\n"
      +"required:\n"
      +"    -i IFILE .................. input molecule file\n"
      +"    -o OFILE .................. output file\n"
      +"\n"
      +"options:\n"
      +"    -imgfmt FMT ............... image format (PNG|JPEG) ["+imgfmt+"]\n"
      +"    -width W .................. width ["+width+"]\n"
      +"    -height H ................. height ["+height+"]\n"
      +"    -kekule ................... kekule\n"
      +"    -arom_gen ................. aromatic general\n"
      +"    -arom_bas ................. aromatic basic\n"
      +"    -showh .................... show Hs\n"
      +"    -transparent .............. transparent (PNG only)\n"
      +"    -maxscale S ............... prevent overscaling small mols ["+maxscale+"]\n"
      +"    -smarts SMARTS ............ pattern match and highlight\n"
      +"    -v[v] ..................... verbose [very]\n"
      +"    -h ........................ this help\n");
    System.exit(1);
  }
  /////////////////////////////////////////////////////////////////////////////
  private static void ParseCommand(String args[])
  {
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-i")) ifile=args[++i];
      else if (args[i].equals("-o")) ofile=args[++i];
      else if (args[i].equals("-width")) width=Integer.parseInt(args[++i]);
      else if (args[i].equals("-height")) height=Integer.parseInt(args[++i]);
      else if (args[i].equals("-imgfmt")) imgfmt=args[++i];
      else if (args[i].equals("-maxscale")) maxscale=Integer.parseInt(args[++i]);
      else if (args[i].equals("-kekule")) kekule=true;
      else if (args[i].equals("-arom_gen")) arom_gen=true;
      else if (args[i].equals("-arom_bas")) arom_bas=true;
      else if (args[i].equals("-showh")) showh=true;
      else if (args[i].equals("-transparent")) transparent=true;
      else if (args[i].equals("-smarts")) smarts=args[++i];
      else if (args[i].equals("-v")) verbose=1;
      else if (args[i].equals("-vv")) verbose=2;
      else if (args[i].equals("-h")) Help("");
      else Help("Unknown option: "+args[i]);
    }
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException
  {
    ParseCommand(args);

    if (ifile==null) Help("Input file required.");

    if (!(new File(ifile).exists())) Help("Non-existent input file: "+ifile);
    MolImporter molReader = new MolImporter(ifile);

    OutputStream ostream=null;
    if (ofile!=null)
      ostream = new FileOutputStream(new File(ofile),false);
    else
      ostream = ((OutputStream)System.out);

    if (verbose>1)
      //System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.getVersion());
      System.err.println("JChem version: "+com.chemaxon.version.VersionInfo.getVersion());

    if (imgfmt.equalsIgnoreCase("jpeg") || imgfmt.equalsIgnoreCase("jpg"))
      imgfmt="jpeg:";
    else
      imgfmt="png:";
    imgfmt+=("h"+height+",w"+width);
    imgfmt+=",maxscale"+maxscale;
    if (kekule)        imgfmt+=",-a";
    else if (arom_gen) imgfmt+=",+a_gen";
    else if (arom_bas) imgfmt+=",+a_bas";
    else               imgfmt+=",-a";
    imgfmt+=(showh?",H_all":",-H");
    if (transparent) imgfmt+=",transbg,#FFFFFF";

    if (verbose>0)
      System.err.println("imgfmt: "+imgfmt);

    Molecule mol=null;
    try { mol=molReader.read(); }
    catch (MolFormatException e)
    {
      System.err.println(e.getMessage());
      System.exit(1);
    }

    if (mol.getDim()==3) mol.clean(2,null,null);

    if (smarts!=null)
    {
      imgfmt+=",setcolors:a1:red:a0:black:b1:red:b0:black";

      MolSearch search = new MolSearch();
      search.setTarget(mol);

      MolHandler smartsReader = null;
      smartsReader = new MolHandler();
      smartsReader.setQueryMode(true);

      try {
        smartsReader.setMolecule(smarts);
      }
      catch (MolFormatException me) {
        me.printStackTrace();
      }

      Molecule query = smartsReader.getMolecule();
      search.setQuery(query);
      boolean ok=false;
      int[][] matchs=null;
      try {
        ok=search.isMatching();
        matchs=search.findAll();
      }
      catch (SearchException se) {
        se.printStackTrace();
      }

      int n_matchs=((matchs!=null)?matchs.length:0);

      System.err.println("smarts matches = "+n_matchs);
      // Atom set 1 represents match atoms, 0 all others:

      Set<Integer> matchbonds_all = new HashSet<Integer>();
      MolBond[] bonds = mol.getBondArray();

      //match atoms will be in set 1, others remain in default set 0
      if (matchs!=null)
      {
        for (int[] match: matchs)
        {
          for (int ia=0; ia<match.length; ++ia)
          {
            if (match[ia]>=0) mol.getAtom(match[ia]).setSetSeq(1); //not impH
          }
        }
        matchbonds_all.addAll(mol2img_utils.getHitBonds(mol, mol2img_utils.getHitAtoms(mol)));

        for (int ib=0;ib<bonds.length;++ib)
          bonds[ib].setSetSeq(matchbonds_all.contains(ib) ? 1:0);
      }
    }

    byte[] data = MolExporter.exportToBinFormat(mol,imgfmt);

    ostream.write(data);
    ostream.close();
    System.exit(0);
  }
  /////////////////////////////////////////////////////////////////////////////
}
