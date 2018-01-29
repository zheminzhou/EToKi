function InteractPlot(div) {
	var self = this;
	this.container = div.css('background', 'white')
		.on('resize', function(e, ui) {
			self.menu.css('width', ui.size.width);
			self.plotView.attr('height', ui.size.height)
				.attr('width', ui.size.width);

			self.plotDim.y0 = ui.size.height - self.margin;
			self.plotDim.x1 = ui.size.width - self.margin;
			var nw = self.plotDim.x1 - self.plotDim.x0;
			self.plotDim.x = nw/self.plotDim.w * (self.plotDim.x - self.plotDim.x0) + self.plotDim.x0;
			self.plotDim.w = nw
			var nh = self.plotDim.y1 - self.plotDim.y0;
			self.plotDim.y = nh/self.plotDim.h * (self.plotDim.y - self.plotDim.y1) + self.plotDim.y1;
			self.plotDim.h = nh;

			self.draw();
		})
	;
	
	this.inSelection = false;
	this.margin = 50;
	this.addMenu();
	this.addTooltip(null, function(d) {
		var data = d[1].map(function(dd) {return dd.id});
		if (data.length > 9) {
			data[9] = '{Other '+(data.length-10)+'}'
			data = data.slice(0, 10)
		}
		return data.join('<br>');
	});
	this.addCanvas();
}

InteractPlot.prototype.addCanvas = function() {
	var self = this;
	this.plotDim = {
		x0 : this.margin, 
		y0 : this.container.height() - this.margin, 
		x1 : this.container.width() - this.margin,
		y1 : this.margin,
		x  : this.margin, 
		y  : this.container.height() - this.margin, 
		w  : this.container.width() - 2*this.margin, 
		h  : 2*this.margin - this.container.height(), 
		dw : 1,
		dh : 1,
	}
	this.container.append('<svg id="plotView"></svg>');
	this.plotView = d3.select(this.container.find('#plotView').selector)
		.attr("height", this.container.height())
		.attr("width", this.container.width())
		.call(d3.drag().on('start', function(d) {
				if (d3.event.sourceEvent.shiftKey) {
					self.dragSelect('start', d3.event.sourceEvent.offsetX, d3.event.sourceEvent.offsetY);
				}
			})
			.on('drag', function(d) {
				if (! self.inSelection) {
					self.plotDim.x = self.plotDim.x+d3.event.dx
					self.plotDim.y = self.plotDim.y+d3.event.dy
					self.draw();
				} else {
					self.dragSelect('drag', d3.event.sourceEvent.offsetX, d3.event.sourceEvent.offsetY);
				}
			})
			.on('end', function(d) {
				self.dragSelect('end');
			})
		)
		.call(d3.zoom().on('zoom', function(d) {
			if (! d3.event.sourceEvent) {
				self.mat.forEach(function(d) {d.style.selected=false;});
				self.plotView.selectAll('.selected').remove();
				return;
			}
			var scale = d3.event.sourceEvent.deltaY < 0 ? 1.1 : 1/1.1;
			if (d3.event.sourceEvent.shiftKey) {
				self.plotDim.dh *= scale;
				self.plotDim.y = d3.event.sourceEvent.offsetY - (d3.event.sourceEvent.offsetY-self.plotDim.y)*scale;
			} else if (d3.event.sourceEvent.ctrlKey) {
				self.plotDim.dh *= scale;
				self.plotDim.y = d3.event.sourceEvent.offsetY - (d3.event.sourceEvent.offsetY-self.plotDim.y)*scale;
				self.plotDim.dw *= scale;
				self.plotDim.x = d3.event.sourceEvent.offsetX - (d3.event.sourceEvent.offsetX-self.plotDim.x)*scale;
			} else {
				self.plotDim.dw *= scale;
				self.plotDim.x = d3.event.sourceEvent.offsetX - (d3.event.sourceEvent.offsetX-self.plotDim.x)*scale;
			}
			self.draw();
		}));
};

InteractPlot.prototype.addTooltip = function(tooltip, func) {
	var self = this;
	func = func ? func : function(d) {return d;};
	self.tooltip = tooltip;
	if (! self.tooltip) {
		self.tooltip = self.container.append('<div id="tooltip"></div>')
					.find('#tooltip')
					.css({
						position : 'fixed',
						width: '100px', 
						height: '300px', 
						opacity: 0.9,
						background : 'cyan', 
						border: '1px solid green', 
					}).hide();
		self.tooltip.on('tooltip:show', function(e, ui) {
			if (ui.data.length == 0) {
				var tt = $(this);
				var dim = [tt.position().left-1, tt.position().top-1, tt.width()+2, tt.height()+2];
				if (!(ui.event.clientX >= dim[0] && ui.event.clientX <= dim[0]+dim[2] && ui.event.clientY >= dim[1] && ui.event.clientY <= dim[1]+dim[3])) {
					$(this).hide();
				}
			} else {
				$(this).show().css({
					left: ui.event.clientX+5,
					top : ui.event.clientY,
				})
				$(this).html(func(ui.data));
			}
		}).on('mouseout', function() {
			$(this).hide();
		})
	}
}


InteractPlot.prototype.addMenu = function() {
	var self = this;
	this.container.addClass('framed')
		.append('<div id="menu" class="tab-handle" style="display:none;position:absolute"></div>')
		.on("contextmenu", function(e, ui) {
			e.preventDefault();
			if (self.menu.css('display') === 'none') {
				self.menu.show(300);
			} else {
				self.menu.hide(300);
			}
		});
	this.menu = this.container.find("#menu")
		.css('width', self.container.width())
		.append('<button id="button-default" class="tab-item-handle" value="default"><span class="glyphicon glyphicon-film"></span></button>')
		.append('<button id="button-resize-vertical" class="tab-item-handle" value="resize-vertical"><span class="glyphicon glyphicon-resize-vertical"></span></button>')
		.append('<button id="button-resize-horizontal" class="tab-item-handle" value="resize-horizontal"><span class="glyphicon glyphicon-resize-horizontal"></span></button>');
	this.menu.find("#button-default").on('click', function(e, ui) {
		self.plotDim.x = self.plotDim.x0;
		self.plotDim.y = self.plotDim.y0;
		self.plotDim.dw = self.plotDim.dh = 1;
		self.draw();
	});
	this.menu.find("#button-resize-vertical").on('click', function(e, ui) {
		self.plotDim.dh *= 1.1;
		self.draw();
	}).on('contextmenu', function(e, ui) {
		e.stopPropagation(); 
		e.preventDefault();
		self.plotDim.dh *= 1/1.1;
		self.draw();
	});
	this.menu.find("#button-resize-horizontal").on('click', function(e, ui) {
		self.plotDim.dw *= 1.1;
		self.draw();
	}).on('contextmenu', function(e, ui) {
		e.stopPropagation(); 
		e.preventDefault();
		self.plotDim.dw *= 1/1.1;
		self.draw();
	})

};


InteractPlot.prototype.dragSelect = function(event, x, y) {
	var self = this;
	if (event == 'end') {
		self.inSelection = false;
		var selectRegion = self.plotView.selectAll('.selection');
		if (selectRegion.size() == 0) {return;}
		var selection = selectRegion.data()[0];
		selectRegion.remove();
		if (selection[0] > selection[2]) {[selection[0], selection[2]] = [selection[2], selection[0]]};
		if (selection[1] > selection[3]) {[selection[1], selection[3]] = [selection[3], selection[1]]};
		var selected = self.plotView.selectAll('.data').filter(function(d) {
			if (d.r === undefined) {
				var dim = [d.x, d.y, d.x+d.w, d.y+d.h];
				var ovlX = Math.min(dim[2], selection[2]) - Math.max(dim[0], selection[0]);
				if (ovlX >= 0) {
					var ovlY = Math.min(dim[3], selection[3]) - Math.max(dim[1], selection[1]);
					if (ovlY >= 0) {
						return true;
					}
				}
			} else {
				if (d.x >= selection[0] && d.x <= selection[2]) {
					if (d.y >= selection[1]-d.r && d.y <= selection[3]+d.r) {
						return true;
					}
				} else if (d.y >= selection[1] && d.y <= selection[3]) {
					if (d.x >= selection[0]-d.r && d.x <= selection[2]+d.r) {
						return true;
					}
				} else {
					var mx = Math.min(Math.abs(selection[0]-d.x), Math.abs(selection[2]-d.x));
					var my = Math.min(Math.abs(selection[1]-d.y), Math.abs(selection[3]-d.y));
					if (mx <= d.r && my <=d.r && d.r*d.r >= mx*mx + my*my ) {
						return true;
					}
				}
			}
			return false;
		});
		if (selected.size() > 0) {
			if (selected.filter(function(d) {return d.style.selected != true;}).size() > 0) {
				selected.each(function(d) {d.style.selected = true;});
			} else {
				selected.each(function(d) {d.style.selected = false;});
			}
			self.draw();
			var datum = []
			self.mat.forEach(function(d) {
				if (d.style.selected === true) {
					datum.push(d.data);
				};
			});
			$(document).trigger("plot:selected", {datum:datum});
		}
	} else if (event == 'start') {
		self.inSelection = true;
		var selection = [[x, y, x, y]];
		self.plotView.selectAll('.selection')
			.data(selection)
			.enter().append('rect')
				.attr('class', 'selection')
				.attr('x', function(d){return d[0];})
				.attr('y', function(d){return d[1];})
				.attr('width', function(d){return d[2]-d[0];})
				.attr('height', function(d){return d[3]-d[1];})
				.attr('fill', 'green')
				.attr('opacity', '0.5');
	} else if (event == 'drag') {
		var selectRegion = self.plotView.selectAll('.selection');
		var selection = selectRegion.data()[0];
		selection[2] = x;
		selection[3] = y;
		if (selection[2] > selection[0]) {
			selectRegion.attr('width', function(d){return d[2]-d[0];});
		} else {
			selectRegion.attr('x', function(d){return d[2];})
				.attr('width', function(d){return d[0]-d[2];});
		}
		if (selection[3] > selection[1]) {
			selectRegion.attr('height', function(d){return d[3]-d[1];});
		} else {
			selectRegion.attr('y', function(d){return d[3];})
				.attr('height', function(d){return d[1]-d[3];});
		}
	}
}
InteractPlot.prototype.prepare = function(data, option) {
	var self = this;
	option = option ? option : {};
	if (option.incremental) {
		self.mat = self.mat.concat(data);
	} else {
		self.mat = data;
	}
	if (option.Xrange) {
		self.Xrange = option.Xrange;
	} else {
		if (data[0].r) {
			self.Xrange = [
				Math.min.apply(Math, data.map(function(d) {return d.x - d.r;})), 
				Math.max.apply(Math, data.map(function(d) {return d.x + d.r;})), 
			];
		} else {
			self.Xrange = [
				Math.min.apply(Math, data.map(function(d) {return d.x;})), 
				Math.max.apply(Math, data.map(function(d) {return d.x + d.w;})), 
			];
		}
	}

	if (option.Yrange) {
		self.Yrange = option.Yrange;
	} else {
		if (data[0].r) {
			self.Yrange = [
				Math.min.apply(Math, data.map(function(d) {return d.y;})), 
				Math.max.apply(Math, data.map(function(d) {return d.y;})), 
			]
		} else {
			self.Yrange = [
				Math.min.apply(Math, data.map(function(d) {return d.y - d.h;})), 
				Math.max.apply(Math, data.map(function(d) {return d.y;})), 
			]
		}
	}
	return self;
}

InteractPlot.prototype.histogram = function(data, option) {
	var self = this;
	option = option ? option : {};
	var x       = option.x        ?       option.x : function(d) {return d[0] ? d[0] : 1;};
	var y       = option.y        ?       option.y : function(d) {return d[1] ? d[1] : 1;};
	var w       = option.w        ?       option.w : function(d) {return d[2] ? d[2] : 1;};
	var h       = option.h        ?       option.h : y;
	var stacked = option.stacked  ? option.stacked : false;
	var style   = option.style    ?   option.style : function(d) {return {stroke:'black', 'stroke-width':0.5, 'fill':(d[3] ? d[3] : 'steelblue')}; };
	
	var mat = data.map(function(d) {return {
						x:x(d), 
						y:y(d), 
						w:w(d), 
						h:h(d),
						style : style(d),
						data : d,
			};
		});
	if (stacked) {
		var stack = {};
		for (var i in mat) {
			var m = mat[i];
			if (stack[m.x]) {
				m.y = stack[m.x] + m.y;
			}
			stack[m.x] = m.y;
			mat.push(m);
		}
	}
	return self.prepare(mat, option);
}
InteractPlot.prototype.heatmap = function(data, option) {
	var self = this;
	var color = option.color ? option.color : ['#f0f0f0', '#000000'];

	var Xmax = Math.max.apply(Math, data.map(function(d) {return d.length;}));
	var Ymax = data.length;
	var Drange = [  Math.min.apply(Math, data.map(function(d) {return Math.min.apply(Math, d);})),
					Math.max.apply(Math, data.map(function(d) {return Math.max.apply(Math, d);})) ];
	var color = d3.scale.linear().domain(Drange)
      .interpolate(d3.interpolateHcl)
      .range(color);

	var mat;
	for (var i in data) {
		var Xinterval = Xmax / data[i].length;
		for (var j in data[i]) {
			var d = data[i][j];
			mat.push({x:Xinterval*j, y:(Ymax-i), w:Xinterval, h:1, data:[i, j, d], style: {fill:color(d)} });
		}
	}
	return self.prepare(mat, option);
}
InteractPlot.prototype.scatter = function(data, option) {
	var self = this;
	option = option ? option : {};
	var x      = option.x        ?      option.x : function(d) {return d[0] ? d[0] : 1;};
	var y      = option.y        ?      option.y : function(d) {return d[1] ? d[1] : 1;};
	var r      = option.r        ?      option.r : function(d) {return d[2] ? d[2] : 1;};
	var style  = option.style    ?  option.style : function(d) {return {stroke:'black', 'stroke-width':0.5, 'fill':'steelblue'}; };
	
	return self.prepare(data.map(function(d) {return {
					x:x(d), 
					y:y(d), 
					r:r(d), 
					style : style(d),
					data : d,
		};
	}), option);
}
InteractPlot.prototype.draw = function(option) {
	var self = this;
	this.plotView.selectAll('*').remove();

	var assign_style = function(item) {
		var style_list = {};
		item.each(function(d) {
			Object.keys(d.style).forEach(function(k) {
					style_list[k] = 1;
				});
		});
		Object.keys(style_list).forEach(function(k) {
			item.attr(k, function(d) {return d.style[k];});
		})
	}

	var Xr = [self.Xrange[0] - 5*(self.Xrange[1] - self.Xrange[0]), self.Xrange[1] + 5*(self.Xrange[1] - self.Xrange[0])];
	var Yr = [self.Yrange[0] - 5*(self.Yrange[1] - self.Yrange[0]), self.Yrange[1] + 5*(self.Yrange[1] - self.Yrange[0])]

	var Xd = [self.plotDim.x, self.plotDim.x+self.plotDim.w*self.plotDim.dw];
	var Yd = [self.plotDim.y, self.plotDim.y+self.plotDim.h*self.plotDim.dh];
	
	Xd = [Xd[0] - 5*(Xd[1] - Xd[0]), Xd[1] + 5*(Xd[1] - Xd[0])];
	Yd = [Yd[0] - 5*(Yd[1] - Yd[0]), Yd[1] + 5*(Yd[1] - Yd[0])];
	
	var x = d3.scaleLinear().domain(Xr).range(Xd);
	var y = d3.scaleLinear().domain(Yr).range(Yd);
	
	var drawing_order = [[-1, 0], [-1, 1], [-1,2]];
	for (var i in self.mat) {
		var m = self.mat[i];
		if (m.r) {
			if (drawing_order[0][0] < 0 ) drawing_order[0][0] = i;
		} else if (m.h) {
			if (drawing_order[1][0] < 0 ) drawing_order[1][0] = i;
		} else {
			if (drawing_order[2][0] < 0 ) drawing_order[2][0] = i;
		}
	}
	drawing_order =	drawing_order.filter(function(d) {return d[0] >= 0})
								.sort(function(d1, d2) {return d1[0] - d2[0]})
								.map(function(d) {return d[1];});
	drawing_order.forEach(function(t) {
		if (t == 0) {
			var item = self.plotView.selectAll('.dot_data').data( self.mat.filter(function(d) {return d.r !== undefined;})
									.map(function(d) {
										return {
											x: x(d.x), 
											y: y(d.y), 
											r: x(d.r) - x(0), 
											style : d.style,
											data : d.data, 
										}})
						).enter().append('g')
							.attr("class", "data");
			var drawing = item.append('circle')
					.attr('cx', function(d) {return d.x;})
					.attr('cy', function(d) {return d.y;})
					.attr('r', function(d) {return d.r;});

			item.filter(function(d) {return d.style.selected === true})
				.append('circle')
				.attr("class", "selected")
				.attr('cx', function(d) {return d.x;})
				.attr('cy', function(d) {return d.y;})
				.attr('r', function(d) {return d.r;})
				.attr("opacity", 0.4)
				.attr("fill", "red")
				.attr("stroke-width", 3)
				.attr("stroke", "black");
			assign_style(drawing);
		} else if (t==1) {
			var item = self.plotView.selectAll('.rect_data').data( self.mat.filter(function(d) {return d.w !== undefined;})
									.map(function(d) {
										return {
											x: x(d.x), 
											y: y(d.y), 
											w: x(d.w) - x(0), 
											h: y(0) - y(d.h), 
											style : d.style,
											data : d.data, 
										}})
						).enter().append('g')
							.attr("class", "data");
			var drawing = item.append('rect')
					.attr("x", function(d) {return d.x;})
					.attr("y", function(d) {return d.y;})
					.attr("width", function(d) {return d.w;})
					.attr("height", function(d) {return d.h;});

			item.filter(function(d) {
				return d.style.selected === true})
				.append('rect')
				.attr("class", "selected")
				.attr("x", function(d) {return d.x;})
				.attr("y", function(d) {return d.y;})
				.attr("width", function(d) {return d.w;})
				.attr("height", function(d) {return d.h;})
				.attr("opacity", 0.4)
				.attr("fill", "red")
				.attr("stroke-width", 3)
				.attr("stroke", "black");
			assign_style(drawing);
		} else {
			item = self.plotView.append("path")
			  .datum( self.mat.filter(function(d) {return (d.r === undefined && d.w === undefined)} ))
			  .attr(d.style)
			  .attr("d", d3.line()
						.x(function(d) {return x(d.x);})
						.y(function(d) {return y(d.y);})
			);
			assign_style(item);
		}
	});

	this.plotView.append('g')
		.attr('transform', "translate(0, "+self.plotDim.y0+")")
		.call(d3.axisBottom(x).ticks(50));
	this.plotView.append('g')
		.attr('transform', "translate(0, "+(self.plotDim.y0+self.plotDim.h)+")")
		.call(d3.axisTop(x).ticks(50));
		
	this.plotView.append('g')
		.attr('transform', "translate("+self.plotDim.x0+",0)")
		.call(d3.axisLeft(y).ticks(50));
	this.plotView.append('g')
		.attr('transform', "translate("+(self.plotDim.x0+self.plotDim.w)+",0)")
		.call(d3.axisRight(y).ticks(50));


	if (self.tooltip) {
		self.container.find('.data').on('mouseover', function(e, ui) {
			self.tooltip.trigger('tooltip:show', {event:e, source:this, data:d3.select(this).data()[0].data});
		})
		.on('mouseout', function(e, ui) {
			self.tooltip.trigger('tooltip:show', {event:e, source:this, data:[]});
		})
	}
	return this;
};
