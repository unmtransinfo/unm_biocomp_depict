package edu.unm.health.biocomp.depict;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.awt.Color;

import org.apache.commons.cli.*; // CommandLine, CommandLineParser, HelpFormatter, OptionBuilder, Options, ParseException, PosixParser

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.util.MolHandler;
import chemaxon.license.LicenseManager;


/**	Utility/test app: Input one molecule, generates PNG or JPEG image.
*/
public class mol2img_app
{
  private static String APPNAME="MOL2IMG";
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

  public static void main(String[] args) throws Exception
  {
    String HELPHEADER =  "MOL2IMG - depict utility";
    Options opts = new Options();
    opts.addOption(Option.builder("i").required().hasArg().desc("Input molecule file").build());
    opts.addOption(Option.builder("o").hasArg().desc("Output file").build());
    opts.addOption(Option.builder("imgfmt").hasArg().desc("Image format (PNG|JPEG) ["+imgfmt+"]").build());
    opts.addOption(Option.builder("width").type(Number.class).hasArg().desc("width ["+width+"]").build());
    opts.addOption(Option.builder("height").type(Number.class).hasArg().desc("height ["+height+"]").build());
    opts.addOption(Option.builder("kekule").desc("Kekule").build());
    opts.addOption(Option.builder("arom_gen").desc("ChemAxon general aromaticity model").build());
    opts.addOption(Option.builder("arom_bas").desc("ChemAxon basic aromaticity model").build());
    opts.addOption(Option.builder("showh").desc("Show hydrogens").build());
    opts.addOption(Option.builder("transparent").desc("transparent (PNG only)").build());
    opts.addOption(Option.builder("maxscale").type(Number.class).hasArg().desc("prevent overscaling small mols ["+maxscale+"]").build());
    opts.addOption(Option.builder("smarts").hasArg().desc("pattern match and highlight").build());
    opts.addOption("v", "verbose", false, "Verbose.");
    opts.addOption("h", "help", false, "Show this help.");
    HelpFormatter helper = new HelpFormatter();
    CommandLineParser clip = new PosixParser();
    CommandLine clic = null;
    try {
      clic = clip.parse(opts, args);
    } catch (ParseException e) {
      helper.printHelp(APPNAME, HELPHEADER, opts, e.getMessage(), true);
      System.exit(0);
    }
    ifile = clic.getOptionValue("i");
    if (clic.hasOption("o")) ofile = clic.getOptionValue("o");
    if (clic.hasOption("imgfmt")) imgfmt = clic.getOptionValue("imgfmt");
    if (clic.hasOption("smarts")) smarts = clic.getOptionValue("smarts");
    if (clic.hasOption("width")) width = (Integer)(clic.getParsedOptionValue("width"));
    if (clic.hasOption("height")) height = (Integer)(clic.getParsedOptionValue("height"));
    if (clic.hasOption("maxscale")) maxscale = (Integer)(clic.getParsedOptionValue("maxscale"));
    if (clic.hasOption("kekule")) kekule = true;
    if (clic.hasOption("arom_gen")) arom_gen = true;
    if (clic.hasOption("arom_bas")) arom_bas = true;
    if (clic.hasOption("showh")) showh = true;
    if (clic.hasOption("transparent")) transparent = true;
    if (clic.hasOption("vv")) verbose = 2;
    else if (clic.hasOption("v")) verbose = 1;
    if (clic.hasOption("h")) {
      helper.printHelp(APPNAME, HELPHEADER, opts, "", true);
      System.exit(0);
    }

    if (!(new File(ifile).exists()))
      helper.printHelp(APPNAME, HELPHEADER, opts, ("Non-existent input file: "+ifile), true);
    MolImporter molReader = new MolImporter(ifile);

    OutputStream ostream=null;
    if (ofile!=null)
      ostream = new FileOutputStream(new File(ofile),false);
    else
      ostream = ((OutputStream)System.out);

    if (verbose>1)
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
      catch (MolFormatException e) {
        e.printStackTrace();
      }

      Molecule query = smartsReader.getMolecule();
      search.setQuery(query);
      boolean ok=false;
      int[][] matchs=null;
      try {
        ok=search.isMatching();
        matchs=search.findAll();
      }
      catch (SearchException e) {
        e.printStackTrace();
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

    byte[] data = MolExporter.exportToBinFormat(mol, imgfmt);

    ostream.write(data);
    ostream.close();
    System.exit(0);
  }
  /////////////////////////////////////////////////////////////////////////////
}
