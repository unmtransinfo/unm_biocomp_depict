package edu.unm.health.biocomp.depict.molalign;

import java.io.*;
import java.util.*;
import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.util.MolHandler;
import chemaxon.sss.search.MolSearch;
import chemaxon.sss.search.SearchException;

/**	Test application for molecular alignment.

	<pre>
	molalign \
	  -i ~/data/2d/nci_quinolines2d_input.sdf -o z.sdf \
	  -smarts '[#7]~1~[#6]~[#6]~[#7]~[#6]~[#6]1'
	</pre>
	@author Jeremy J Yang
*/
public class molalign
{
  private static void Help(String msg)
  {
    if (!msg.equals("")) System.err.println(msg);
    System.err.println(
      "usage: molalign\n"+
      "  required:\n"+
      "          -i <in_mol_file>\n"+
      "          -o <out_mol_file> ... w/ aligned 2D\n"+
      "\n"+
      "          -smarts <smarts>\n"+
      "       or\n"+
      "          -mces ... align by MCES to 1st mol (BROKEN BY API CHANGE)\n"+
      "  options:\n"+
      "          -ifmt <fmt_spec>\n"+
      "          -ofmt <fmt_spec>\n"+
      "          -v    ... verbose\n"
    );
    System.exit(1);
  }

  /////////////////////////////////////////////////////////////////////////////
  public static void main(String[] args)
    throws IOException,SearchException
  {
    if (args.length==0) Help("");
    String ifile="";
    String smarts="";
    String ofile="";
    String ifmt="";
    String ofmt=null;
    boolean mces=false;
    int verbose=0;
    for (int i=0;i<args.length;++i)
    {
      if (args[i].equals("-i")) { ifile=args[++i]; }
      else if (args[i].equals("-smarts")) { smarts=args[++i]; }
      else if (args[i].equals("-mces")) { mces=true; }
      else if (args[i].equals("-o")) { ofile=args[++i]; }
      else if (args[i].equals("-ifmt")) { ifmt=args[++i]; }
      else if (args[i].equals("-ofmt")) { ofmt=args[++i]; }
      else if (args[i].equals("-v")) { verbose=1; }
      else {
        Help("Bad option: "+args[i]);
      }
    }
    if (ifile.equals("") || ofile.equals(""))
      Help("-i and -o required");
    if (smarts.equals("") && !mces)
      Help("-mces or -smarts required");

    MolImporter molReader;
    if (ifile.equals("-"))
    {
      if (ifmt.equals("")) Help("-ifmt required with \"-i -\"");
      molReader=new MolImporter(System.in,ifmt);
    }
    else
    {
      molReader=new MolImporter(ifile);
    }

    MolExporter molWriter=null;
    if (ofile.equals("-"))
    {
      if (ofmt==null) Help("-ofmt required with \"-o -\"");
      molWriter=new MolExporter(System.out,ofmt);
    }
    else
    {
      if (ofmt==null)
        ofmt=MFileFormatUtil.getMostLikelyMolFormat(ofile);
      molWriter=new MolExporter(new FileOutputStream(ofile),ofmt);
    }

    if (verbose>0)
      //System.err.println("JChem version: "+chemaxon.jchem.version.VersionInfo.getVersion());
      System.err.println("JChem version: "+com.chemaxon.version.VersionInfo.getVersion());

    ArrayList<Molecule> mols = new ArrayList<Molecule>();
    Molecule mol;
    while ((mol=molReader.read())!=null)
    {
      mol.aromatize(MoleculeGraph.AROM_GENERAL);
      if (mol.getDim()!=2) mol.clean(2,null,null);
      mols.add(mol);
    }
    if (mols.size()==0)
    {
      System.err.println("ERROR: No molecules!");
      System.exit(1);
    }

    if (!smarts.equals(""))
    {
      MolHandler smartsReader = new MolHandler();
      smartsReader.setQueryMode(true);
      try { smartsReader.setMolecule(smarts); }
      catch (MolFormatException e) {
        System.err.println(e.getMessage());
        System.exit(1);
      }
      MolSearch search = new MolSearch();
      Molecule qmol = smartsReader.getMolecule();
      search.setQuery(qmol);

      search.setTarget(mols.get(0));
      int [] hitatoms = search.findFirst();
      if (hitatoms==null)
      {
        System.err.println("ERROR: no match on first molecule.  Quitting.");
        System.exit(1);
      }
    }

    int n_aligned=0;
    if (mces)
    {
      //n_aligned=molalign_utils.AlignToMCES(mols);
      System.err.println("ERROR: MCES alignment broken due to API change...");
      System.exit(1);
    }
    else
    {
      n_aligned=molalign_utils.AlignToSmarts(mols,smarts);
    }

    for (int i=0;i<mols.size();++i)
    {
      if (verbose>0)
        System.err.println(mols.get(i).getName());
      molWriter.write(mols.get(i));
    }
    System.err.println("mols processed: "+mols.size()+" mols aligned: "+n_aligned);
    System.exit(0);
  }
}
