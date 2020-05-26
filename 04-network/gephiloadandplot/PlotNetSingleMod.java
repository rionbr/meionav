/*
author: Rion Brattig Correia
description:
  loads a graphml file and plots it using the gephi toolkit.
  largely based from https://github.com/klarman-cell-observatory/forceatlas2
*/
import com.itextpdf.text.PageSize;

import org.gephi.graph.api.Graph;
import org.gephi.graph.api.GraphController;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.Node;
import org.gephi.graph.api.Table;
import org.gephi.graph.api.Column;
import org.gephi.appearance.api.AppearanceController;
import org.gephi.appearance.api.AppearanceModel;
import org.gephi.io.importer.api.Container;
import org.gephi.io.importer.api.EdgeDirectionDefault;
import org.gephi.io.importer.api.ImportController;
import org.gephi.io.processor.plugin.DefaultProcessor;
import org.gephi.io.exporter.api.ExportController;
import org.gephi.io.exporter.preview.PDFExporter;
import org.gephi.io.exporter.preview.PNGExporter;
import org.gephi.preview.api.PreviewController;
import org.gephi.preview.api.PreviewModel;
import org.gephi.preview.api.PreviewProperty;
import org.gephi.preview.types.DependantColor;
import org.gephi.preview.types.EdgeColor;
import org.gephi.layout.spi.Layout;
import org.gephi.layout.plugin.rotate.RotateLayout;
import org.gephi.project.api.ProjectController;
import org.gephi.project.api.Workspace;
import org.openide.util.Lookup;
//
import java.awt.Color;
import java.io.*;
import java.util.*;
import com.itextpdf.text.PageSize;


public class PlotNetSingleMod {

    private static Map<String, Arg> argsMap = new LinkedHashMap<>();

    private static void addArg(String flag, String description, boolean not_boolean, Object defaultValue) {
        argsMap.put("--" + flag.toLowerCase(), new Arg(flag, description, not_boolean, "" + defaultValue));
    }

    private static void addArg(String flag, String description, boolean not_boolean) {
        argsMap.put("--" + flag.toLowerCase(), new Arg(flag, description, not_boolean, null));
    }

    private static String getArg(String flag) {
        Arg a = argsMap.get("--" + flag.toLowerCase());
        return a != null ? a.value : null;
    }

    private static class Arg {
        String flag;
        String description;
        boolean not_boolean;
        String defaultValue;
        String value;

        private Arg(String flag, String description, boolean not_boolean, String defaultValue) {
            this.flag = flag;
            this.description = description;
            this.not_boolean = not_boolean;
            this.defaultValue = defaultValue;
            if (defaultValue != null) {
                this.value = defaultValue;
            }
        }

        public String toString() {
            return flag;
        }
    }

    public static void main(String[] args) throws IOException {
        
        addArg("input", "Input graph in one of Gephi input file formats https://gephi.org/users/supported-graph-formats/", true);
        addArg("coords", "Input file with node x,y coordinates", true);
        addArg("output", "Output .pdf file.", true);
        addArg("directed", "Whether input graph is undirected.", false, false);
        addArg("rotate", "Rotate the graph 'x' degrees.", true);
        addArg("highlight", "Which attribute to highlight.", true);

        for (int i = 0; i < args.length; i++) {
            Arg a = argsMap.get(args[i].toLowerCase());
            if (a == null) {
                System.err.println("Unknown argument " + args[i]);
                System.exit(1);
            }
            String value = a.not_boolean ? args[++i] : "true";
            a.value = value;
        }

        File coordsFile = null;

        // Input File
        File file = new File(getArg("input"));
        if (!file.exists()) {
            System.err.println(file + " not found.");
            System.exit(1);
        }

        // Coordinates
        if (getArg("coords") != null) {
            coordsFile = new File(getArg("coords"));
            if (!coordsFile.exists()) {
                System.err.println(coordsFile + " not found.");
                System.exit(1);
            }
        }
        

        // Layout Rotate Angle
        Double rotate = null;
        if (getArg("rotate") != null) {
            rotate = Double.parseDouble(getArg("rotate"));
        }

        // Highlight attribute
        String highlight = null;
        if (getArg("highlight") != null) {
            highlight = getArg("highlight");
        }

        //Init a project - and therefore a workspace
        ProjectController pc = Lookup.getDefault().lookup(ProjectController.class);
        pc.newProject();
        Workspace workspace = pc.getCurrentWorkspace();

        //Get models and controllers for this new workspace - will be useful later
        GraphController graphController = Lookup.getDefault().lookup(GraphController.class);
        GraphModel graphModel = graphController.getGraphModel();
        ImportController importController = Lookup.getDefault().lookup(ImportController.class);

        //Import file
        Container container = importController.importFile(file);
        Graph graph;
        if (!getArg("directed").equalsIgnoreCase("true")) {
            container.getLoader().setEdgeDefault(EdgeDirectionDefault.UNDIRECTED);
            graph = graphModel.getUndirectedGraph();
        } else {
            container.getLoader().setEdgeDefault(EdgeDirectionDefault.DIRECTED);
            graph = graphModel.getDirectedGraph();
        }
        //See if graph is well imported
        //System.out.println("Nodes: " + g.getNodeCount());
        //System.out.println("Edges: " + g.getEdgeCount());

        //Append imported data to GraphAPI
        importController.process(container, new DefaultProcessor(), workspace);

        // Node Positioning from coords file
        if (coordsFile != null) {
            Map<Object, Node> idToNode = new HashMap<>();
            for (Node n : graph.getNodes()) {
                idToNode.put(n.getId(), n);

            }
            BufferedReader br = new BufferedReader(new FileReader(coordsFile));
            String sep = "\t";
            String s = br.readLine();
            for (String test : new String[]{"\t", ","}) {
                if (s.indexOf(test) != -1) {
                    sep = test;
                    break;
                }
            }
            List<String> header = Arrays.asList(s.split(sep));
            int idIndex = header.indexOf("id");
            int xIndex = header.indexOf("x");
            int yIndex = header.indexOf("y");
            int zIndex = header.indexOf("z");
            while ((s = br.readLine()) != null) {
                String[] tokens = s.split(sep);
                String id = tokens[idIndex];
                Node n = idToNode.get(id);
                if (n != null) {
                    n.setX(Float.parseFloat(tokens[xIndex]));
                    n.setY(Float.parseFloat(tokens[yIndex]));
                } else {
                    System.err.println(id + " not found");
                }
            }
            br.close();
        }


        //Model
        PreviewModel previewModel = Lookup.getDefault().lookup(PreviewController.class).getModel();


        // Attributes
        /*
        System.out.println("--Node Attributes--");
        for (Column col : graphModel.getNodeTable()) {
            System.out.println("- " + col);
        }
        */

        // Rotate layout
        if (rotate != null) {
            RotateLayout rotateLayout = new RotateLayout(null, rotate);
            rotateLayout.setGraphModel(graphModel);
            //
            rotateLayout.initAlgo();
            if (rotateLayout.canAlgo()) {
                rotateLayout.goAlgo();
            }
            rotateLayout.endAlgo();
        }

        // Node
        previewModel.getProperties().putValue(PreviewProperty.NODE_BORDER_WIDTH, 0.0f);
        previewModel.getProperties().putValue(PreviewProperty.NODE_BORDER_COLOR, new DependantColor(Color.decode("#000000")));
        previewModel.getProperties().putValue(PreviewProperty.SHOW_NODE_LABELS, Boolean.FALSE);

        for(Node node : graphModel.getGraph().getNodes()) {
            //System.out.println(node.getId() + " : " + node.getLabel() ); // + " : " + node.getAttribute("new-DM-phenotype"));
            //System.out.println(node.getAttribute(Column ""))

            // Highlight node based on column
            if (!"none".equals(highlight)) {
                if (node.getAttribute(highlight) == Boolean.TRUE) {
                    node.setColor(Color.decode("#d62728")); //red
                } else {
                    node.setColor(Color.decode("#000000")); //blue
                }
            }
            // Node Size
            node.setSize(6.0f);
        }

        //Edge
        previewModel.getProperties().putValue(PreviewProperty.EDGE_COLOR, new EdgeColor(EdgeColor.Mode.MIXED));
        previewModel.getProperties().putValue(PreviewProperty.EDGE_OPACITY, new Float(5.0f));
        previewModel.getProperties().putValue(PreviewProperty.EDGE_THICKNESS, new Float(0.1f));

        //Export PDF
        ExportController ec = Lookup.getDefault().lookup(ExportController.class);
        //
        PDFExporter pdfExporter = (PDFExporter) ec.getExporter("pdf");
        pdfExporter.setPageSize(PageSize.A4);
        pdfExporter.setWorkspace(workspace);
        pdfExporter.setLandscape(Boolean.TRUE);
        ec.exportFile(new File(getArg("output")), pdfExporter);
    }
}
