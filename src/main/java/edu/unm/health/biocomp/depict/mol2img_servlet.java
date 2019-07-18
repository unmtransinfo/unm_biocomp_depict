package edu.unm.health.biocomp.depict;

import java.io.*;
import java.text.*;
import java.util.*;
import java.util.regex.*;
import java.awt.Color;
import javax.servlet.*;
import javax.servlet.http.*;

import chemaxon.formats.*;
import chemaxon.struc.*;
import chemaxon.sss.search.*;
import chemaxon.util.MolHandler;
import chemaxon.license.LicenseManager;

import edu.unm.health.biocomp.util.http.*;

/**	Generates PNG or JPEG image for inline display.
	Input smiles or MDL molfile.
	Use existing 2D unless specified otherwise.
	Note that molfiles and other connection table formats are expected to be gzipped and base64-encoded.
	<UL>
	<LI> Depicts smarts 
	<LI> Smarts matching 
	<LI> Atom colors control
	</UL>
	to do: [ ] superatoms (e.g. OEt, CF3)
	to do: [ ] reduce size of OH2 etc.g.
	Configuration note: use CATALINA_OPTS="-Djava.awt.headless=true" for PNG/JPEG java.awt calls to work in
	Tomcat environment.
	<br>
	@author Jeremy J Yang
*/
public class mol2img_servlet extends HttpServlet
{
  /////////////////////////////////////////////////////////////////////////////
  public void doPost(HttpServletRequest request,HttpServletResponse response)
	throws IOException,ServletException,MolFormatException
  {
    ResourceBundle rb=ResourceBundle.getBundle("LocalStrings",request.getLocale());
    LicenseManager.refresh();
    String imgfmt=request.getParameter("imgfmt");
    if (imgfmt!=null && (imgfmt.equalsIgnoreCase("jpeg") || imgfmt.equalsIgnoreCase("jpg")))
    {
      imgfmt="jpeg:";
      response.setContentType("image/jpeg");
    }
    else
    {
      imgfmt="png:";
      response.setContentType("image/png");
    }
    String mdl=request.getParameter("mdlcode");
    if (mdl==null) mdl=request.getParameter("mdl");
    String mrv=request.getParameter("mrvcode");
    if (mrv==null) mrv=request.getParameter("mrv");
    String smiles=request.getParameter("smiles");
    if (smiles==null) smiles=request.getParameter("smi");
    if (smiles==null) smiles=request.getParameter("smicode");
    boolean smi_is_cx=(request.getParameter("smi_is_cx")!=null);
    String fcode=request.getParameter("fcode");
    String format=request.getParameter("format");
    if (format!=null && (format.equalsIgnoreCase("sd") || format.equalsIgnoreCase("mol") || format.equalsIgnoreCase("mdl"))) format="sdf";
    if (format!=null && format.equalsIgnoreCase("sdf") && fcode!=null) mdl=fcode;

    String width=request.getParameter("w");
    if (width==null || width.length()==0) width="300";
    String height=request.getParameter("h");
    if (height==null || height.length()==0) height="200";
    imgfmt+=("h"+height+",w"+width);

    String maxscale=request.getParameter("maxscale");	//Max magnification, prevent overscaling small molecules.
    if (maxscale==null) maxscale="28"; //default
    imgfmt+=",maxscale"+maxscale;

    String kekule=request.getParameter("kekule");
    if (kekule!=null && kekule.equalsIgnoreCase("true"))
      imgfmt+=",-a";
    String arom_gen=request.getParameter("arom_gen");
    if (arom_gen!=null && arom_gen.equalsIgnoreCase("true"))
      imgfmt+=",+a_gen";
    String arom_bas=request.getParameter("arom_bas");
    if (arom_bas!=null && arom_bas.equalsIgnoreCase("true"))
      imgfmt+=",+a_bas";
    String showh=request.getParameter("showh");
    if (showh!=null && showh.equalsIgnoreCase("true"))
      imgfmt+=",H_all";
    else
      imgfmt+=",-H";
    String hideh=request.getParameter("hideh");
    if (hideh!=null && hideh.equalsIgnoreCase("true"))
      imgfmt+=",H_off";
    String colormode=request.getParameter("mode");
    String bgcolor=request.getParameter("bgcolor");
    if (bgcolor!=null)
      imgfmt+=","+bgcolor;  // e.g. png:w100,#FFFF00
    else if (colormode!=null && colormode.equals("bow"))
      imgfmt+=",mono";
    else if (colormode!=null && colormode.equals("cob"))
      imgfmt+=",#000000";
    else
      imgfmt+=",#FFFFFF";
    String lonepairs=request.getParameter("lonepairs");
    if (lonepairs!=null && lonepairs.equalsIgnoreCase("true"))
      imgfmt+=",lp";
    String showmaps=request.getParameter("showmaps");
    if (showmaps!=null && showmaps.equalsIgnoreCase("true"))
      imgfmt+=",amap";
    String transparent=request.getParameter("transparent");
    if (transparent!=null && transparent.equalsIgnoreCase("true"))
      imgfmt+=",transbg,#FFFFFF";
    String smarts=request.getParameter("smarts");
    if (smarts==null) smarts=request.getParameter("smartscode");
    String smartses=request.getParameter("smartses");
    boolean smilesmatch=(request.getParameter("smilesmatch")!=null);  //Meaning: smarts is a smiles, so relax syntax, ignore Hs, e.g. "[nH]".
 
    String matchmrv_code=request.getParameter("matchmrv");  // alternate match/highlight method
    Molecule qmol=null;
    if (matchmrv_code!=null)
    {
      qmol=MolImporter.importMol(matchmrv_code,"base64:gzip:mrv");
    }

    // atomcolors is list of color codes [0-9]+.
    String atomcolors=request.getParameter("atomcolors");

    Molecule mol=null;
    if (mdl!=null) 
    {
      //byte[] mdlgz=Base64Decoder.decodeToBytes(mdl) ; //Now JChem can do this...
      try { mol=MolImporter.importMol(mdl,"base64:gzip:mol"); }
      catch (MolFormatException e) { throw new MolFormatException(e.getMessage()+"\n"+mdl,e); }
    }
    else if (mrv!=null) 
    {
      //byte[] mrvgz=Base64Decoder.decodeToBytes(mrv) ; //Now JChem can do this...
      try { mol=MolImporter.importMol(mrv,"base64:gzip:mrv"); }
      catch (MolFormatException e) { throw new MolFormatException(e.getMessage()+"\n"+mrv,e); }
    }
    else if (smiles!=null) 
    {
      try { mol=MolImporter.importMol(smiles,(smi_is_cx?"cxsmiles:":"smiles:")); }
      catch (MolFormatException e) { throw new MolFormatException(e.getMessage()+"\n"+smiles,e); }
    }
    else if (format!=null && fcode!=null)
    {
      try { mol=MolImporter.importMol(fcode,"base64:gzip:"+format); }
      catch (MolFormatException e) { throw new MolFormatException(e.getMessage()+"\n"+fcode,e); }
    }
    if (mol==null) throw new IOException("ERROR: mol==null.  smiles: \""+smiles+"\"");

    ArrayList<String> smartslist = new ArrayList<String>();
    if (smartses!=null)
    {
      for (String s: Pattern.compile("\\t").split(smartses))
        smartslist.add(s);
    }
    else if (smarts!=null)
    {
      smartslist.add(smarts);
    }

    if (smartslist.size()>0 && (new MolSearch()).isLicensed())
    {
      try { mol2img_utils.ColorBySmarts(mol,smartslist,smilesmatch); } catch (SearchException e) { } //depict anyway
      imgfmt+=",setcolors:a0:black:b0:black:b1:red:a1:red:b2:red"; //required after JChem5_2_0:
    }
    else if (qmol!=null)
    {
      try { mol2img_utils.ColorByQMol(mol,qmol); } catch (SearchException e) { } //depict anyway
      imgfmt+=",setcolors:a0:black:b0:black:b1:red:a1:red:b2:red"; //required after JChem5_2_0:
    }

    String clearqprops=request.getParameter("clearqprops");
    if (clearqprops!=null && clearqprops.equalsIgnoreCase("true"))
    {
      for (MolAtom atom:mol.getAtomArray()) atom.clearQProps();
    }

    // According to forum topic 5126, the valence display (e.g. "(v0)")
    // can be turned off in upcoming version 5.4, but not till then.
    // Atomcolors idxs correspond with MolAtom::setSetSeq() and mrvSetSeq vals.
    if (atomcolors!=null)
    {
      ArrayList<String> colors = new ArrayList<String>(); // format RRGGBB
      for (int i=0;i<10;++i)
      {
        String color=request.getParameter("color"+i);
        if (color!=null && color.charAt(0)=='#') color=color.substring(1);
        colors.add(color);
      }
      mol2img_utils.ColorByNumbers(mol,atomcolors,colors);
      imgfmt+=",setcolors"; //required after JChem5_2_0:
      //imgfmt+=":a0:000000:a1:440044:a2:660066:a3:660066:a4:880088:a5:880088:a6:CC0044:a7:CC0044:a8:FF0000:a9:FF0000"; // default
      for (int i=0;i<colors.size();++i)
      {
        if (colors.get(i)==null) continue;
        else if (colors.get(i).length()!=6)
        {
          //System.err.println("DEBUG: bad color="+colors.get(i));
          continue; //ERROR
        }
        imgfmt+=":a"+i+":"+colors.get(i);
      }
    }
    if (mol.getDim()==3) mol.clean(2,null,null);

    // output:
    //byte[] data=mol.toBinFormat(imgfmt); //old way
    byte[] data;
    try { data = MolExporter.exportToBinFormat(mol,imgfmt); }
    catch (IOException e) { throw new IOException(e.getMessage()+"\n"+("DEBUG: mol==null: "+((mol==null)?"YES":"NO")),e); }
    ServletOutputStream ostream=response.getOutputStream();
    ostream.write(data);
    ostream.close();
  }
  /////////////////////////////////////////////////////////////////////////////
  public void doGet(HttpServletRequest request,HttpServletResponse response)
	throws IOException,ServletException,MolFormatException
  {
    doPost(request,response);
  }
}
