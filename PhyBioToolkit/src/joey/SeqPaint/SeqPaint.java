/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package joey.SeqPaint;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Container;
import java.awt.Font;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.util.HashMap;
import java.util.Map;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JScrollBar;
import javax.swing.JSplitPane;
import javax.swing.WindowConstants;
import org.biojava.bio.BioRuntimeException;
import org.biojava.bio.gui.sequence.AlignmentRenderer;
import org.biojava.bio.gui.sequence.GUITools;
import org.biojava.bio.gui.sequence.LabelledSequenceRenderer;
import org.biojava.bio.gui.sequence.MultiLineRenderer;
import org.biojava.bio.gui.sequence.RulerRenderer;
import org.biojava.bio.gui.sequence.SequenceRenderContext;
import org.biojava.bio.gui.sequence.SymbolSequenceRenderer;
import org.biojava.bio.gui.sequence.TranslatedSequencePanel;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.SimpleAlignment;
import org.biojava.bio.symbol.SymbolList;
import org.biojava3.core.sequence.ProteinSequence;

/**
 *
 * @author zhongxf
 */
public class SeqPaint 
{
    TranslatedSequencePanel tsp = new TranslatedSequencePanel();
    
    //LabelledSequenceRenderer for each AlignmentRenderer
    LabelledSequenceRenderer labRen1, labRen2, labRen3, labRen;
   
    //AlignmentRenderer to hold each sequence
    AlignmentRenderer render1, render2, render3, render;
    //MultiLineRenderer to allow display of multiple tracks in the TranslatedSequencePanel
    MultiLineRenderer multi = new MultiLineRenderer();
    
    //SymbolSequenceRenderer to handle display of the sequence symbols - only one instance is needed
    SymbolSequenceRenderer symbol = new SymbolSequenceRenderer();
    //RulerRenderer to display sequence coordinates
    RulerRenderer ruler = new RulerRenderer();
    //The width in pixels of the of the label in the LabelledSequenceRenderer 
    int labelWidth = 50;
    //The height in pixels of the of the label in the LabelledSequenceRenderer 
    int labelHeight = 15;
    JScrollBar scrollBar;
    
    public void paintSeq(Sequence[] nrSeqs) 
    {
        SymbolSequenceRenderer symbol = new SymbolSequenceRenderer() 
        {
            Paint outline = null;
            double depth = 15;

            public void paint(Graphics2D g2, SequenceRenderContext context) 
            {
                Rectangle2D prevClip = g2.getClipBounds();
                AffineTransform prevTransform = g2.getTransform();
                g2.setPaint(outline);
                Font font = context.getFont();
                Rectangle2D maxCharBounds = font.getMaxCharBounds(g2.getFontRenderContext());
                double scale = context.getScale();
                if (scale >= (maxCharBounds.getWidth() * 0.3) && scale >= (maxCharBounds.getHeight() * 0.3)) 
                {
                    double xFontOffset = 0.0;
                    double yFontOffset = 0.0;

                    xFontOffset = maxCharBounds.getCenterX() * 0.25;
                    yFontOffset = -maxCharBounds.getCenterY() + (depth * 0.5);

                    SymbolList seq = context.getSymbols();
                    SymbolTokenization toke = null;
                    try 
                    {
                        toke = seq.getAlphabet().getTokenization("token");
                    } catch (Exception ex) 
                    {
                        throw new BioRuntimeException(ex);
                    }
                    Location visible = GUITools.getVisibleRange(context, g2);
                    for (int sPos = visible.getMin(); sPos <= visible.getMax(); sPos++) 
                    {
                        double gPos = context.sequenceToGraphics(sPos);
                        String s = "?";
                        try 
                        {
                            s = toke.tokenizeSymbol(seq.symbolAt(sPos));
                        } 
                        catch (Exception ex) 
                        {
                            // We'll ignore the case of not being able to tokenize it
                        }

                        //Start of the modifications -------------------------------
                        //Make sure the text is uppercase
                        s = s.toUpperCase();
                        //Set the color according to the nucleotide for the background
                        if (s.equals("A"))      {g2.setColor(new Color(255, 140, 105));} 
                        else if (s.equals("T")) {g2.setColor(new Color(238, 238, 0));} 
                        else if (s.equals("G")) {g2.setColor(new Color(176, 226, 255));} 
                        else if (s.equals("C")) {g2.setColor(new Color(151, 251, 152));} 
                        else                    {g2.setColor(new Color(230, 230, 250)); }

                        //Plot the rectangle to frame the nucleotide symbol
                        g2.fill(new Rectangle2D.Double((gPos + xFontOffset) - 1.5, 0, tsp.getScale(), labelHeight));
                        //Set the colour for the text
                        g2.setColor(new Color(83, 83, 83));
                        //End of the modifications ---------------------------------

                        g2.drawString(s, (float) (gPos + xFontOffset), (float) yFontOffset);
                    }
                }
                g2.setTransform(prevTransform);
            }
        };

        //Use the Map to create a new SimpleAlignment
        Map<String, Sequence> list = new HashMap();
        for (int i =0; i < nrSeqs.length; i ++)
        {
            list.put(nrSeqs[i].getName(), nrSeqs[i]);
        }
//        list.put(nrSeqs[0].getName(), nrSeqs[0]);
//        list.put(nrSeqs[1].getName(), nrSeqs[1]);
//        list.put(nrSeqs[2].getName(), nrSeqs[2]);
        SimpleAlignment ali = new SimpleAlignment((Map) list);
         multi.addRenderer(ruler);

         //for (int j =0; j <nrSeqs.length; j++) ===display alignment from bottom to top
        for (int j =nrSeqs.length-1; j >=0; j--) //===display alignment from top to bottom
        {
            //Instantiate the AlignmentRenderer
        render = new AlignmentRenderer();
        //Set the label for the AlignmentRenderer
        render.setLabel(ali.getLabels().get(j));
        //Set the renderer for the AlignmentRenderer
        render.setRenderer(symbol);
        //Instantiate the LabelledSequenceRenderer
        labRen = new LabelledSequenceRenderer(labelWidth, labelHeight);
        //Set the name of the sequence as the label in the LabelledSequenceRenderer
        labRen.addLabelString(render.getLabel().toString());
        //Put the AlignmentRenderer in the LabelledSequenceRenderer
        labRen.setRenderer(render);
        //Add the alignment renderers to the MultiLineRenderer
        multi.addRenderer(labRen);
       
        //Add the ruler to the MultiLineRenderer
//        multi.addRenderer(ruler);
        }
         multi.addRenderer(ruler);


        //Set the sequence in the TranslatedSequencePanel
        tsp.setSequence((SymbolList) ali);
        //Set the background colour of the TranslatedSequencePanel
        tsp.setOpaque(true);
        tsp.setBackground(Color.white);
        //Set the renderer for the TranslatedSequencePanel
        tsp.setRenderer(multi);

        //Create a scrollbar and add an adjustment listener
        scrollBar = new JScrollBar(JScrollBar.HORIZONTAL, 0, 0, 0, 100);
        scrollBar.addAdjustmentListener(
                new AdjustmentListener() {
                    public void adjustmentValueChanged(AdjustmentEvent e) {
                        //Get the absolute position of the scroll bar
                        double scrollBarValue = e.getValue();
                        //Get the position of the scroll bar relative to the maximum value
                        double scrollBarRatio = scrollBarValue / scrollBar.getMaximum();
                        //Calculate the new position of the first base to be displayed
                        double pos = scrollBarRatio * (tsp.getSequence().length() - ((tsp.getWidth() - labelWidth) / tsp.getScale()));
                        //Set the new SymbolTranslation for the TranslatedSequencePanel
                        tsp.setSymbolTranslation((int) Math.round(pos));
                    }
                });

        // Set up the display
        JFrame jf = new JFrame("NRseqPaint");
        Container con = jf.getContentPane();
        con.setLayout(new BorderLayout());
        con.add(tsp, BorderLayout.CENTER);
        con.add(scrollBar, BorderLayout.SOUTH);
        jf.setSize(400, 170);
        jf.setLocation(100, 100);
        jf.setVisible(true);
        jf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

    }

    public void paintSeq(ProteinSequence[] aaSeqs) throws IllegalSymbolException 
    {
        SymbolSequenceRenderer symbol = new SymbolSequenceRenderer() 
        {
            Paint outline = null;
            double depth = 20;
            public void paint(Graphics2D g2, SequenceRenderContext context) 
            {
                Rectangle2D prevClip = g2.getClipBounds();
                AffineTransform prevTransform = g2.getTransform();
                
                g2.setPaint(outline);
                Font font = context.getFont();
                Rectangle2D maxCharBounds = font.getMaxCharBounds(g2.getFontRenderContext());
                double scale = context.getScale();
                if (scale >= (maxCharBounds.getWidth() * 0.3) && scale >= (maxCharBounds.getHeight() * 0.3)) 
                {
                    double xFontOffset = 0.0;
                    double yFontOffset = 0.0;

                    xFontOffset = maxCharBounds.getCenterX() * 0.25;
                    yFontOffset = -maxCharBounds.getCenterY() + (depth * 0.5);

                    SymbolList seq = context.getSymbols();
                    SymbolTokenization toke = null;

                    try 
                    {
                        toke = seq.getAlphabet().getTokenization("token");
                    } catch (Exception ex) 
                    {
//                  throw new BioRuntimeException(ex);
                    }
                    Location visible = GUITools.getVisibleRange(context, g2);

                    for (int sPos = visible.getMin(); sPos <= visible.getMax(); sPos++) 
                    {
                        double gPos = context.sequenceToGraphics(sPos);
                        String s = "?";
                        try 
                        {
                            s = toke.tokenizeSymbol(seq.symbolAt(sPos));
                        } catch (Exception ex) 
                        {
                            // We'll ignore the case of not being able to tokenize it
                        }

                        //Start of the modifications -------------------------------
                        //Make sure the text is uppercase
                        s = s.toUpperCase();
                        System.out.println(s);

                        //Set the color according to the nucleotide for the background
                        if (s.equals("A")) {
                            g2.setColor(Color.darkGray);
                        } else if (s.equals("V")) {
                            g2.setColor(Color.darkGray);
                        } else if (s.equals("L")) {
                            g2.setColor(Color.darkGray);
                        } else if (s.equals("I")) {
                            g2.setColor(Color.darkGray);
                        } else if (s.equals("M")) {
                            g2.setColor(Color.darkGray);
                        } else if (s.equals("S")) {
                            g2.setColor(Color.green);
                        } else if (s.equals("T")) {
                            g2.setColor(Color.green);
                        } else if (s.equals("N")) {
                            g2.setColor(Color.green);
                        } else if (s.equals("Q")) {
                            g2.setColor(Color.green);
                        } else if (s.equals("D")) {
                            g2.setColor(Color.red);
                        } else if (s.equals("E")) {
                            g2.setColor(Color.red);
                        } else if (s.equals("H")) {
                            g2.setColor(Color.blue);
                        } else if (s.equals("R")) {
                            g2.setColor(Color.blue);
                        } else if (s.equals("K")) {
                            g2.setColor(Color.blue);
                        } else if (s.equals("F")) {
                            g2.setColor(Color.magenta);
                        } else if (s.equals("Y")) {
                            g2.setColor(Color.magenta);
                        } else if (s.equals("W")) {
                            g2.setColor(Color.magenta);
                        } else if (s.equals("P")) {
                            g2.setColor(Color.cyan);
                        } else if (s.equals("G")) {
                            g2.setColor(Color.cyan);
                        } else {
                            g2.setColor(Color.black);
                        }

                        g2.fill(new Rectangle2D.Double((gPos + xFontOffset) - 1.5, 0, tsp.getScale(), labelHeight));

//            g2.fill(new Rectangle2D.Double((gPos + xFontOffset)-1.5, 0, tsp.getScale(), labelHeight ));
                        //Set the colour for the text
                        g2.setColor(new Color(83, 83, 83));
                        //End of the modifications ---------------------------------

                        g2.drawString(s, (float) (gPos + xFontOffset), (float) yFontOffset);
                    }

                }
                g2.setTransform(prevTransform);
            }
        };

        Map<String, SymbolList> list = new HashMap();
        for (int l = 0; l < aaSeqs.length; l++) {
            list.put(aaSeqs[l].getAccession().toString(), ProteinTools.createProtein(aaSeqs[l].toString()));
        }


        SimpleAlignment ali = new SimpleAlignment((Map) list);
        String[] titleArray = new String[aaSeqs.length];
        //Add the alignment renderers to the MultiLineRenderer
        for (int ilen = 0; ilen < aaSeqs.length; ilen++) {
           
            //Instantiate the AlignmentRenderer
            render = new AlignmentRenderer();
            //Set the label for the AlignmentRenderer
            render.setLabel(ali.getLabels().get(ilen));
            //Set the renderer for the AlignmentRenderer
            render.setRenderer(symbol);
            


            //Instantiate the LabelledSequenceRenderer
            labRen = new LabelledSequenceRenderer(labelWidth, labelHeight);
            //Set the name of the sequence as the label in the LabelledSequenceRenderer
            labRen.addLabelString(render.getLabel().toString());
            //Put the AlignmentRenderer in the LabelledSequenceRenderer
            labRen.setRenderer(render);
            

             titleArray[ilen] = render.getLabel().toString();
//        //## if you want to render sequence with a label, use this
        multi.addRenderer(labRen);  
            // render sequences without labels
//            multi.addRenderer(render);
           
        }

        multi.addRenderer(ruler);

        //Set the sequence in the TranslatedSequencePanel
        tsp.setSequence((SymbolList) ali);
        //Set the background colour of the TranslatedSequencePanel
        tsp.setOpaque(true);
        tsp.setBackground(Color.white);
        //Set the renderer for the TranslatedSequencePanel
        tsp.setRenderer(multi);
        
       
        // Set up the display
        JFrame jf = new JFrame("aaSeqPaint");
        jf.setTitle("hello world");
        Container con = jf.getContentPane();
        con.setLayout(new BorderLayout());
        
        JList titleList = new JList(titleArray);
//        titleList.setSize(80, 400);
//        titleList.setFixedCellHeight(20);
        
        JSplitPane jSplitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT);
        jSplitPane.setLeftComponent(titleList);
        jSplitPane.setRightComponent(tsp);
//        con.add(tsp, BorderLayout.CENTER);
        con.add(jSplitPane, BorderLayout.CENTER);
        JScrollBar jsc = new JScrollBar();
        //con.add(jsc, BorderLayout.SOUTH);
        jf.setSize(400, 400);
        jf.setLocation(200, 200);
        jf.setVisible(true);
        jf.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

    }

    public TranslatedSequencePanel getSeqPanel() 
    {
        return tsp;
    }

    
}
