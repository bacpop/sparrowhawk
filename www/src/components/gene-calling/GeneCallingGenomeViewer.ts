import React, { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react';
import { createViewState, JBrowseApp } from '@jbrowse/react-app2';
import makeWorkerInstance from '@jbrowse/react-app2/esm/makeWorkerInstance';
import MgnifyGeneViewerPlugin from 'mgnify-jbrowse/dist/components/GeneViewer/jbrowse/plugin';
import { FeaturePanel } from 'mgnify-jbrowse/dist/components/GeneViewer/components/FeaturePanel';
import { fetchFirstFaiRef, queryGffRegion } from 'mgnify-jbrowse/dist/components/GeneViewer/gff';
import { useGeneViewerSessionSync } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useGeneViewerSessionSync';
import { useGeneViewerTrackRefresh } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useGeneViewerTrackRefresh';
import { useGeneViewerHideDrawer } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useGeneViewerHideDrawer';
import { useGeneViewerResizeSync } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useGeneViewerResizeSync';
import { useGeneViewerClickHandler } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useGeneViewerClickHandler';
import { useGeneViewerSelection } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useGeneViewerSelection';
import { useJBrowseVisibleRegion } from 'mgnify-jbrowse/dist/components/GeneViewer/hooks/useJBrowseVisibleRegion';
import type { GffFeature } from 'mgnify-jbrowse/dist/components/GeneViewer/gff';
import type PluginManager from '@jbrowse/core/PluginManager';

interface GeneCallingGenomeViewerProps {
    fileName: string;
    fastaUrl: string;
    faiUrl: string;
    gziUrl: string;
    gffUrl: string;
    csiUrl: string;
    heightPx?: number;
}

const AMR_COLOR = '#dc2626';
const GENE_COLOR = '#6b7280';
const SELECTED_COLOR = '#2563eb';

function getFeatureValue(feature: any, key: string): string {
    if (!feature) return '';
    const lower = key.toLowerCase();
    if (typeof feature.get === 'function') {
        const direct = feature.get(key) ?? feature.get(lower);
        if (direct != null && direct !== '') return String(direct);
        const attrs = feature.get('attributes');
        if (attrs && typeof attrs === 'object') {
            const attr = attrs[key] ?? attrs[lower];
            if (attr != null && attr !== '') return String(attr);
        }
    }
    const data = feature.data ?? feature;
    if (data && typeof data === 'object') {
        const direct = data[key] ?? data[lower];
        if (direct != null && direct !== '') return String(direct);
        const attrs = data.attributes;
        if (attrs && typeof attrs === 'object') {
            const attr = attrs[key] ?? attrs[lower];
            if (attr != null && attr !== '') return String(attr);
        }
    }
    return '';
}

function featureId(feature: any): string {
    return getFeatureValue(feature, 'ID') || getFeatureValue(feature, 'id') || '';
}

class GeneCallingJBrowsePlugin extends MgnifyGeneViewerPlugin {
    name = 'GeneCallingJBrowsePlugin';

    install(pluginManager: PluginManager): void {
        super.install(pluginManager);
        pluginManager.jexl.addFunction('getGeneCallingColor', (feature: any) => {
            const selected = typeof window !== 'undefined' ? String((window as any).selectedGeneId ?? '').trim() : '';
            const id = featureId(feature).trim();
            if (selected && id === selected) return SELECTED_COLOR;
            return getFeatureValue(feature, 'amr_unit_label') ? AMR_COLOR : GENE_COLOR;
        });
    }
}

function buildAssemblyConfig(props: GeneCallingGenomeViewerProps) {
    return {
        name: props.fileName,
        displayName: props.fileName,
        sequence: {
            type: 'ReferenceSequenceTrack',
            trackId: 'ReferenceSequenceTrack',
            adapter: {
                type: 'BgzipFastaAdapter',
                fastaLocation: { uri: props.fastaUrl },
                faiLocation: { uri: props.faiUrl },
                gziLocation: { uri: props.gziUrl },
            },
        },
    };
}

function buildGeneTrackConfig(props: GeneCallingGenomeViewerProps) {
    const labelJexl = ["Name", "gene", "locus_tag", "ID"].map(f => `get(feature,'${f}')`).join(' || ');
    return {
        type: 'FeatureTrack',
        trackId: 'gene_features',
        name: 'Genes',
        assemblyNames: [props.fileName],
        category: ['Annotations'],
        adapter: {
            type: 'Gff3TabixWithEssentialityAdapter',
            gffGzLocation: { uri: props.gffUrl },
            index: {
                indexType: 'CSI',
                location: { uri: props.csiUrl },
            },
            featureJoinAttribute: 'ID',
        },
        displays: [
            {
                displayId: 'gene_features-LinearBasicDisplay',
                id: 'gene_features-LinearBasicDisplay',
                type: 'LinearBasicDisplay',
                height: 280,
                renderer: {
                    type: 'SvgFeatureRenderer',
                    color1: 'jexl:getGeneCallingColor(feature)',
                    color2: 'jexl:getGeneCallingColor(feature)',
                    labels: { name: 'jexl:' + labelJexl },
                },
            },
        ],
        visible: true,
    };
}

function buildSessionConfig(fileName: string, refName: string, end: number, geneTrack: any) {
    return {
        name: 'Gene Calling session',
        widgets: { BaseFeatureWidget: { type: 'BaseFeatureWidget', disabled: true } },
        views: [
            {
                type: 'LinearGenomeView',
                configuration: { header: { disable: true, hidden: true } },
                displayedRegions: [{ refName, start: 0, end, assemblyName: fileName }],
                tracks: [
                    {
                        id: geneTrack.trackId,
                        type: 'FeatureTrack',
                        configuration: geneTrack.trackId,
                        minimized: false,
                        visible: true,
                        displays: geneTrack.displays,
                    },
                ],
            },
        ],
    };
}

function patchFeaturePanelHeading(root: HTMLElement | null): void {
    if (!root) return;
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    let node = walker.nextNode();
    while (node) {
        if (node.textContent === 'Feature Details') {
            node.textContent = 'Feature details';
        }
        node = walker.nextNode();
    }
}

function patchFeatureContextMenus(viewState: any): void {
    try {
        const views = viewState?.session?.views ?? [];
        for (const view of views) {
            for (const track of view.tracks ?? []) {
                for (const display of track.displays ?? []) {
                    if (typeof display.contextMenuItems !== 'function' || (display as any).__amrMenuPatched) continue;
                    const original = display.contextMenuItems.bind(display);
                    Object.defineProperty(display, 'contextMenuItems', {
                        configurable: true,
                        value: () => original().filter((item: any) => item?.label !== 'Open feature details'),
                    });
                    (display as any).__amrMenuPatched = true;
                }
            }
        }
    } catch {
        // If the JBrowse model blocks patching, the side panel still provides working details.
    }
}

export function GeneCallingGenomeViewer(props: GeneCallingGenomeViewerProps) {
    const [viewState, setViewState] = useState<any>(null);
    const [error, setError] = useState<string | null>(null);
    const [selectedGeneId, setSelectedGeneId] = useState<string | null>(null);
    const [genesInView, setGenesInView] = useState<GffFeature[]>([]);
    const genesInViewRef = useRef<GffFeature[]>([]);
    genesInViewRef.current = genesInView;
    const lastTableSelectionTimeRef = useRef(0);
    const containerRef = useRef<HTMLDivElement | null>(null);
    const sidePanelRef = useRef<HTMLDivElement | null>(null);
    const heightPx = props.heightPx ?? 600;

    const assemblyConfig = useMemo(() => buildAssemblyConfig(props), [props.fileName, props.fastaUrl, props.faiUrl, props.gziUrl]);
    const geneTrackConfig = useMemo(() => buildGeneTrackConfig(props), [props.fileName, props.gffUrl, props.csiUrl]);

    useEffect(() => {
        let cancelled = false;
        async function init() {
            try {
                setError(null);
                const first = await fetchFirstFaiRef(props.faiUrl);
                const sessionConfig = buildSessionConfig(props.fileName, first.refName, first.length, geneTrackConfig);
                const state = createViewState({
                    config: {
                        assemblies: [assemblyConfig],
                        tracks: [{ ...geneTrackConfig, visible: true }],
                        defaultSession: { ...sessionConfig, name: 'defaultSession' },
                    },
                    plugins: [GeneCallingJBrowsePlugin],
                    makeWorkerInstance,
                });
                try {
                    const session = (state as any).session;
                    session.showWidget = () => undefined;
                    session.addWidget = () => undefined;
                } catch {
                    // ignore
                }
                patchFeatureContextMenus(state);
                if (!cancelled) setViewState(state);
            } catch (e: any) {
                if (!cancelled) setError(e?.message ?? String(e));
            }
        }
        init();
        return () => { cancelled = true; };
    }, [assemblyConfig, geneTrackConfig, props.fileName, props.faiUrl]);

    useEffect(() => {
        if (viewState) patchFeatureContextMenus(viewState);
    }, [viewState]);

    useLayoutEffect(() => {
        if (typeof window !== 'undefined') {
            (window as any).selectedGeneId = selectedGeneId ?? '';
        }
    }, [selectedGeneId]);

    const visibleRegion = useJBrowseVisibleRegion(viewState, 500);

    useEffect(() => {
        if (!visibleRegion) return;
        let cancelled = false;
        const { refName, start, end } = visibleRegion;
        const regionLen = end - start;
        const buffer = Math.max(0, Math.floor(regionLen * 0.25));
        window.setTimeout(async () => {
            try {
                const features = await queryGffRegion({
                    gffUrl: props.gffUrl,
                    csiUrl: props.csiUrl,
                    refName,
                    start: Math.max(0, start - buffer),
                    end: end + buffer,
                    featureTypes: ['CDS'],
                });
                if (!cancelled) setGenesInView(features);
            } catch (e: any) {
                if (!cancelled) setError(e?.message ?? String(e));
            }
        }, 150);
        return () => { cancelled = true; };
    }, [visibleRegion?.refName, visibleRegion?.start, visibleRegion?.end, props.gffUrl, props.csiUrl]);

    const resolveToLocusTag = useCallback((id: string, features: GffFeature[]) => {
        const norm = String(id).trim();
        const match = features.find(f => String(f.attributes?.ID ?? f.id ?? '').trim() === norm);
        return String(match?.attributes?.ID ?? match?.id ?? norm);
    }, []);

    const emptyEssentialityIndex = useMemo(() => new Map<string, string>(), []);
    const { selectedFeatures } = useGeneViewerSelection(selectedGeneId, genesInView, 'ID', false, emptyEssentialityIndex, undefined);

    useGeneViewerSessionSync({ viewState, lastTableSelectionTimeRef, setSelectedGeneId, joinAttr: 'ID' });
    useGeneViewerTrackRefresh(viewState, selectedGeneId, selectedGeneId, emptyEssentialityIndex, false);
    useGeneViewerClickHandler({
        viewState,
        containerRef,
        genesInViewRef,
        lastTableSelectionTimeRef,
        setSelectedGeneId,
        resolveToLocusTag,
        joinAttribute: 'ID',
    });
    useGeneViewerHideDrawer(viewState, containerRef);
    useGeneViewerResizeSync(viewState, containerRef);

    useLayoutEffect(() => {
        patchFeaturePanelHeading(sidePanelRef.current);
    }, [selectedFeatures]);

    return React.createElement('div', { style: { width: '100%', border: '1px solid #d1d5db', borderRadius: 8, overflow: 'hidden' } },
        error ? React.createElement('div', { style: { padding: 12, background: '#fef2f2', borderBottom: '1px solid #fecaca', color: '#991b1b' } }, error) : null,
        React.createElement('div', { style: { display: 'grid', gridTemplateColumns: '1fr 320px', width: '100%' } },
            React.createElement('div', { ref: containerRef, style: { width: '100%', minWidth: 0, minHeight: heightPx, maxHeight: heightPx, overflow: 'hidden' } },
                viewState ? React.createElement(JBrowseApp, { viewState }) : React.createElement('div', { style: { padding: 12, color: '#6b7280' } }, 'Loading JBrowse...')
            ),
            React.createElement('div', { ref: sidePanelRef, style: { borderLeft: '1px solid #d1d5db', height: heightPx, overflowY: 'auto', overflowX: 'hidden' } },
                React.createElement(FeaturePanel, { features: selectedFeatures, essentiality: null })
            )
        )
    );
}
