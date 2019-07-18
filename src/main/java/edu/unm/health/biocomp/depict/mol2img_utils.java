package edu.unm.health.biocomp.depict;

import java.io.*;
import java.text.*;
import java.util.*;
import java.util.regex.*;
import java.awt.Color;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.*;
import chemaxon.sss.search.*;
import chemaxon.util.MolHandler;

/**	Molecule image utilities, static methods.
	<br>
	@author Jeremy J Yang
*/
public class mol2img_utils
{
  /////////////////////////////////////////////////////////////////////////////
  /**	Set sequence numbers for atoms matched by smarts for subsequent
	coloring.
  */
  public static void ColorBySmarts(Molecule mol,String smarts,
	boolean smilesmatch)
	throws MolFormatException,SearchException
  {
    ArrayList<String> smartslist = new ArrayList<String>();
    smartslist.add(smarts);
    ColorBySmarts(mol,smartslist,smilesmatch);
  }
  /////////////////////////////////////////////////////////////////////////////
  public static void ColorBySmarts(Molecule mol,String smarts)
	throws MolFormatException,SearchException
  {
    ColorBySmarts(mol,smarts,false);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Set sequence numbers for atoms matched by smarts for subsequent
	coloring.
  */
  public static void ColorBySmarts(Molecule mol,ArrayList<String> smartslist)
	throws MolFormatException,SearchException
  {
    ColorBySmarts(mol,smartslist,false);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Set sequence numbers for atoms matched by smarts for subsequent
	coloring.
	<br>
	Example for smilesmatch: <br>
	mol smiles = Cc1c(\N=C2/C=CC(=O)c3ccccc23)c(=O)n(-c4ccccc4)n1C <br>
	match smiles = O=c1cc[nH]n1c1ccccc1 <br>
	<br>
	@param	smilesmatch	regard smarts as smiles; relax syntax, ignore H's (e.g. [nH])
  */
  public static void ColorBySmarts(Molecule mol,ArrayList<String> smartslist,
	boolean smilesmatch)
	throws MolFormatException,SearchException
  {
    // Must aromatize here so smarts work correctly,
    mol.dearomatize();  // May help for some, maybe "C1CC=c2[nH]c3ccccc3cc2C1"
    mol.aromatize(MoleculeGraph.AROM_GENERAL);

    MolSearch search = new MolSearch();
    //MolHandler molHand = new MolHandler();

    //Default set is 0.
    //mol.setAtomSetSeqs(0); // Atom set 1 match, 0 non-match.
    //mol.setBondSetSeqs(0); // Bond set 2 match, 0 non-match.
    // The API behavior is a bit weird...
    // For some reason using bond set 1 does not work?

    HashSet<Integer> matchbonds_all=new HashSet<Integer>(); //all smarts,matches
    for (String sma:smartslist)
    {
      search.setTarget(mol);
      Molecule qmol;
      if (smilesmatch)
      {
        qmol=MolImporter.importMol(sma,"smiles:");
        // Remove exp-H's from query molecule, ignore imp-H's for permissive match.
        qmol.implicitizeHydrogens(MolAtom.ALL_H);
        qmol.dearomatize();  // May help for some, maybe "C1CC=c2[nH]c3ccccc3cc2C1"
        qmol.aromatize(MoleculeGraph.AROM_GENERAL);
        MolSearchOptions searchOpts = new MolSearchOptions(SearchConstants.SUBSTRUCTURE);
        searchOpts.setImplicitHMatching(SearchConstants.IMPLICIT_H_MATCHING_DISABLED);
        search.setSearchOptions(searchOpts);
      }
      else
      {
        qmol=MolImporter.importMol(sma,"smarts:");
        //molHand.setMolecule(sma);
        //molHand.setQueryMode(true);
        //qmol=molHand.getMolecule();
      }
      search.setQuery(qmol);
      boolean ok=false;
      int[][] matchs=null;
      ok=search.isMatching();
      matchs=search.findAll(); //may throw SearchException

      if (matchs!=null)
      {
        //Each match is an array of atom indices.
        for (int[] match: matchs)
        {
          for (int ia=0; ia<match.length; ++ia)
            if (match[ia]>=0) mol.getAtom(match[ia]).setSetSeq(1); //not impH

// Discontinued!  MolHandler.getNonHitBonds()
//          HashSet<Integer> nhbonds_this = new HashSet<Integer>();
//          for (Object bond: MolHandler.getNonHitBonds(qmol,mol,match))
//            nhbonds_this.add(mol.indexOf((MolBond)bond));
//          for (int ib=0;ib<mol.getBondCount();++ib)
//            if (!matchbonds_all.contains(ib) && !nhbonds_this.contains(ib)) matchbonds_all.add(ib);
        }
      }
    }
    matchbonds_all.addAll(getHitBonds(mol, getHitAtoms(mol)));
    for (int ib:matchbonds_all)
      mol.getBond(ib).setSetSeq(matchbonds_all.contains(ib) ? 2:0);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Set sequence numbers for atoms matched by query molecule for subsequent
	coloring.

  */
  public static void ColorByQMol(Molecule mol,Molecule qmol)
	throws SearchException
  {
    // Must aromatize here to match correctly,
    mol.dearomatize();  // May help for some, maybe "C1CC=c2[nH]c3ccccc3cc2C1"
    mol.aromatize(MoleculeGraph.AROM_GENERAL);
    MolSearch search=new MolSearch();
    search.setQuery(qmol);
    search.setTarget(mol);
    boolean ok=false;
    int[][] matchs=null;
    ok=search.isMatching();
    matchs=search.findAll(); //may throw SearchException
    if (matchs==null) return;

    //Each match is an array of atom indices.
    HashSet<Integer> matchbonds_all=new HashSet<Integer>(); //all matches
    for (int[] match: matchs)
    {
      for (int ia=0; ia<match.length; ++ia)
        if (match[ia]>=0) mol.getAtom(match[ia]).setSetSeq(1); //not impH

// Discontinued!  MolHandler.getNonHitBonds()
//      HashSet<Integer> nhbonds_this = new HashSet<Integer>();
//      for (Object bond: MolHandler.getNonHitBonds(qmol,mol,match))
//        nhbonds_this.add(mol.indexOf((MolBond)bond));
//      for (int ib=0;ib<mol.getBondCount();++ib)
//        if (!matchbonds_all.contains(ib) && !nhbonds_this.contains(ib)) matchbonds_all.add(ib);
    }

    matchbonds_all.addAll(getHitBonds(mol, getHitAtoms(mol)));

    for (int ib:matchbonds_all)
      mol.getBond(ib).setSetSeq(matchbonds_all.contains(ib) ? 2:0);
  }
  /////////////////////////////////////////////////////////////////////////////
  /**	Set sequence numbers denoting atom colors.
  	@param	acolorstr	color codes for each atom in order
  	@param	colors	color palette, hex encoded
  */
  public static boolean ColorByNumbers(Molecule mol,String acolorstr,ArrayList<String> colors)
  {
    int natoms = mol.getAtomCount();
    char [] tmp = acolorstr.toCharArray();
    if (tmp.length!=natoms)
    {
      //System.err.println("DEBUG: tmp.length="+tmp.length+" != natoms="+natoms);
      return false;     //ERROR
    }
    Integer[] acolors = new Integer[natoms];
    for (int i=0;i<natoms;++i)
    {
      int j = Integer.parseInt(""+tmp[i]);
      if (j<0 || j>9 || colors.get(j)==null)
      {
        return false;   //ERROR
      }
      acolors[i]=j;
    }
    Integer [] colors_rgb = new Integer [10];
    for (int i=0;i<10;++i)
    {
      if (colors.get(i)==null) colors_rgb[i]=null;
      else
      {
        colors_rgb[i]=Integer.parseInt(colors.get(i),16);   //parse as hex
        if (colors_rgb[i]>0xFFFFFF) return false;       //ERROR
      }
    }
    // Define atom sets.
    MolAtom[] atoms = mol.getAtomArray();
    for (int i=0;i<natoms;++i)
    {
//      if (atoms[i].hasAromaticBond())
//      {
//        for (int j=0;j<atoms[i].getBondCount();++j)
//        { atoms[i].getBond(j).setSetSeq(acolors[i]); }
//      }
      atoms[i].setSetSeq(acolors[i]);
      //System.err.println("DEBUG: atom("+i+") color="+acolors[i]);
    }
    //Color bonds (seems to be needed for aromatic circles.)
    //Using bond coloring to clear unwanted half bond colors.
//    MolBond[] bonds = mol.getBondArray();
//    for (MolBond bond:bonds)
//    {
//      MolAtom a1=bond.getAtom1();
//      MolAtom a2=bond.getAtom2();
//      if (a1.getSetSeq()==a2.getSetSeq()) bond.setSetSeq(a1.getSetSeq());
//    }

    return true;
  }

  /** From ChemAxon forum post 61743. */
  public static Set<Integer> getHitAtoms(Molecule mol)
  {
    Set<Integer> hitAtoms = new HashSet<Integer>();
    for (int i=0; i<mol.getAtomCount(); i++)
    {
      if (mol.getAtom(i).getSetSeq()>0)
      {
        hitAtoms.add(i);
      }
    }
    return hitAtoms;
  }
  public static Set<Integer> getHitBonds(Molecule mol, Set<Integer> hitAtoms)
  {
    Set<Integer> hitBonds = new HashSet<Integer>();
    for (int i=0; i<mol.getBondCount(); i++)
    {
      MolBond bond = mol.getBond(i);
      if (bond.getAtom1().getSetSeq()>0 && bond.getAtom2().getSetSeq()>0)
      {
        hitBonds.add(i);
      }
    }
    return hitBonds;
  }
}
