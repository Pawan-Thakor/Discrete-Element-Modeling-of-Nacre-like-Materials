%Plot Displacement for the Trusses

function Plot_Displacement(p)
    us=1:2:2*p.Num_Nodes-1; 
    vs=2:2:2*p.Num_Nodes;
    XX=p.Displacement(us);
    YY=p.Displacement(vs); 
    dispNorm=max(sqrt(XX.^2+YY.^2)); 
    scaleFact=2*dispNorm;
    
    h = figure;
    plot(1)

    figure(1); 
    hold on;
    plot(p.Node_coordinate(:,1),p.Node_coordinate(:,2),'k.')
    hold on;
    axis equal;  
        for e=1:p.Num_Elements       
            el_node = p.Element_Nodes(e, 1:2);
            node_xy = p.Node_coordinate(el_node, :);
            plot(node_xy(:,1), node_xy(:,2),'k--');
        end
     node_new = p.Node_coordinate + scaleFact*reshape(p.Displacement,2,p.Num_Nodes)';
     plot(node_new(:,1), node_new(:,2),'o','MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',10);
     hold on;
     axis equal;
         for e=1:p.Num_Elements       
             el_node = p.Element_Nodes(e, 1:2);
             node_xy = node_new(el_node, :);
             plot(node_xy(:,1), node_xy(:,2),'k-','LineWidth',2);
         end
        
end